#!/usr/bin/env python3
# coding=utf-8

import sys, os, re
from argparse import ArgumentParser, FileType
import numpy as np
import powerlaw

"""
Estimate power law distribution params (alpha and xmin) based on (a set of) stored
error simulation histogram(s). The result is an approximation of the actual fitted
power law params from the simulation since the stored historgram(s) have only the
number of observations in each log-likelihood ratio bin, not the actual log-
likelihood ratios for each error observation. Using the default bin step size of
0.5, this should still produce a reasonably close approximation of the actual
log-likelihood ratio distribution, just somewhat less granular.
"""

__version__ = "1.0.1"


def parse_histogram_file(log_file):
    """
    Parse the error histogram data from the given file.
    """    
    with open(log_file, 'r') as infile:
        # Need the EOF position
        infile.seek(0, os.SEEK_END)
        eof = infile.tell()
        infile.seek(0)

        # Get ready to parse the file
        found = False
        block_start_pos = block_end_pos = last_block_start_pos = last_block_end_pos = infile.tell()
        prev_block_line_count = block_line_count = 0
        while not found:
            pos = infile.tell()
            line = infile.readline().strip('\n')
            #sys.stderr.write("%s\n" % line)
            if re.match('iteration', line):
                # Found the start of a histogram block
                block_start_pos = pos
                block_line_count = 0
            elif pos == eof:
                # We've reached the end of the file
                if block_line_count >= prev_block_line_count:
                    # We have a complete block
                    last_block_start_pos = block_start_pos
                    last_block_end_pos = pos
                else:
                    last_block_end_pos = block_end_pos
                # Stop looking when we've parsed the whole file
                found = True
            elif re.fullmatch('', line):
                # Complete histogram blocks are separated by an empty line
                #sys.stderr.write("%s\n" % (pos))
                block_end_pos = pos
                last_block_start_pos = block_start_pos
                prev_block_line_count = block_line_count
            elif re.match("\d+", line):
                # Histogram line. Keep a count.
                block_line_count += 1
                
            else:
                # Something unexpected. Do nothing.
                pass
    
        #sys.stderr.write("%s\t%s\n" % (last_block_start_pos, last_block_end_pos))

        # Parse the last complete block for histogram bins and counts
        _hist = []
        _bins = []
        bin_width = 0
        i = 0
        infile.seek(last_block_start_pos)        
        while infile.tell() < last_block_end_pos:
            line = infile.readline().strip('\n')
            if re.match("\d+", line):
                # log-likelihood bin lines start with a digit
                #sys.stderr.write("%s\n" % (line))
                if i == 1:
                    # Calculate the bin width on the second bin
                    bin_width = float(line[0]) - _bins[0]
                line = line.split('\t')
                # Store the bin and observed number of errors in this bin
                _hist.append(int(line[2]))
                _bins.append(float(line[0]))
                i += 1

        hist = np.array(_hist)
        bins = np.array(_bins)
        #sys.stderr.write("%s\n%s\t%s\n" % (hist, bins, bin_width))
    return hist, bins, bin_width


def combine_hists_and_bins(hist_all, bins_all, hist_i, bins_i, bin_width):
    """
    Combine two sets of error counts and bins.
    """
    if len(bins_i) > len(bins_all):
        bins_all = np.append(bins_all, np.arange(max(bins_all),
                                                 max(bins_i)+bin_width,
                                                 bin_width)
        )
        hist_all = sum_unequal_hists(hist_all, hist_i)
    elif len(bins_all) > len(bins_i):
        hist_all = sum_unequal_hists(hist_i, hist_all)
    else:
        hist_all = hist_all + hist_i

    return hist_all, bins_all


def sum_unequal_hists(shorter_arr, longer_arr):
    """
    Return the sum of two numpy histograms with different maximum bins.
    """
    return longer_arr + np.append(shorter_arr, np.zeros(len(longer_arr) - len(shorter_arr)))
    

def convert_hist_to_likelihood_ratios(hist, bins):
    """
    Create a vector of log-likelihood ratios from the histogram of
    error counts in each bin.
    """
    # Need the overal number of error observations
    rlen = np.sum(hist)
    # Initialize the return vector as an np.array with default values
    # so we can address by index (avoids needing to use a nested loop)
    ret = np.arange(rlen, dtype=float)
    pos = 0
    for i in range(0, len(bins)):
        # Loop over bins and store N values equal to the bin value
        # where N is the count stored in the histogram
        ret[pos:pos+hist[i]] = bins[i]
        pos += hist[i]
    #sys.stderr.write("%s\n" % (ret))
    return ret.tolist()


def main():
    parser = ArgumentParser(description="")
    parser.add_argument('input_files', type=str, nargs='+',
                        help='Per-iteration error stats file(s) containing stats to parse as input for powerlaw parameter fitting.')
    parser.add_argument('--output_file', '-o',  type=FileType('w'), nargs='?',
                        default='simulation_error_model_output.txt',
                        help='Output file name/path. Default=simulation_error_model_output.txt')
    parser.add_argument('--stats_file', '-s',  type=FileType('w'), nargs='?',
			default='simulation_error_stats.tsv',
                        help='File into which error stats from each file will be summarized. Default=simulation_error_stats.tsv')    
    
    args = parser.parse_args()

    # Read input file(s) and combine observation counts in each bin, if multiple files are given.
    sys.stderr.write("\n################ Reading Input Files ################\n")
    hist_all = np.array([])
    bins_all = np.array([])
    i = 0
    for infile in args.input_files:
        # Find the last complete set of histogram bins and extract the counts for each bin into hist_err.
        hist, bins, bin_width = parse_histogram_file(infile)
        # Roll these data into the total data
        if (i > 0):
            hist_all, bins_all = combine_hists_and_bins(hist_all, bins_all, hist, bins, bin_width)
        else:
            hist_all, bins_all = (hist, bins)
        #sys.stderr.write("%s\n%s\n" % (hist_all, bins_all))
        i += 1

    # Extrapolate a vector of log-likelihood ratios from the histogram of counts
    error_lr_list = convert_hist_to_likelihood_ratios(hist_all, bins_all)
        
    # Fit the power law distribution
    sys.stderr.write("\n################ Fitting Error Model ################\n")
    pl_fit = powerlaw.Fit(error_lr_list, method="KS")
    
    # Apply a linear xmin correction based on sample-size. 3.585 is the intercept term derived
    # from fitting a linear model to the difference between actual and estimated xmin values
    # over many sets of values sampled from powerlaw distributions with known parameters. We have
    # shown this correction to reduce the bias in xmin estimates across a wide range of sample sizes.
    corrected_xmin = ( pl_fit.xmin * (3.585 / np.log(len(error_lr_list))) ) + 3

    # Write results to output files
    args.output_file.write("Total errors in simulations: %d\nMin log-likelihood-ratio for errors: %.3f\nMax log-likelihood ratio for errors: %.3f\nPowerlaw fitting method: %s\nNumber of observations used in powerlaw fitting: %d\nFitted powerlaw alpha: %.5f\nFitted powerlaw xmin: %.5f\nCorrected powerlaw xmin: %.5f\nFitted powerlaw minimum xmin distance (D): %.3f\nFitted powerlaw sigma: %.5f\n" % (len(error_lr_list), min(error_lr_list), max(error_lr_list), pl_fit.fit_method, pl_fit.n, pl_fit.alpha, pl_fit.xmin, corrected_xmin, pl_fit.D, pl_fit.sigma) )

        
if __name__ == '__main__':
    main()


# coding=utf-8

import sys, os, re
#from argparse import ArgumentParser, Namespace, _ArgumentGroup, _SubParsersAction
import argparse
from LRphase.SimulatePhasedData import *
from LRphase.InputData import *
from LRphase import urls
from pysam import AlignmentFile, VariantFile
import numpy as np
from scipy.stats import nbinom
import powerlaw
import time
from LRphase.PhaseSet import powlaw_modified, _estimate_prior_probabilities

__version__ = "1.1.1"

# TO-DO:
## The main phasing function should be a wrapper that loops over samples
# in the input VCF. Content of the current function should be moved out
# and called once per sample in the VCF, unless args.one_sample is given.
#@profile
def phasing(args):
    """Main function for phasing mode

    Args:
        args:
    """
    sys.stderr.write('\nAnalyzing inputs using phasing mode...\n\n')
    #sys.stderr.write("%s\n" % args)
    # Input data initialization currently ignores much of the content of args. May be a good idea to fix this.
    # TO-DO: Need to actually handle VCF samples properly instead of relying on args.one_sample.
    input_data = InputData(args.output_directory_name,
                           sample = args.one_sample,
                           ignore_samples = args.ignore_samples,
                           reference_sequence_input = args.reference_genome,
                           reference_sequence_input_assembly = args.reference_assembly,
                           max_threads = args.threads,
                           quiet = args.quiet_mode,
                           silent = args.silent_mode
    )
    input_data.add_haplotype_information(args.vcf_file_name)
    input_data.add_reference_sequence(args.reference_genome,
                                      reference_sequence_input_assembly = args.reference_assembly,
                                      output_directory = args.output_directory_name)
    try:
        input_data.add_reads(args.long_read_inputs,
                             sample = args.one_sample,
                             quiet = args.quiet_mode,
                             silent = args.silent_mode,
                             threads = args.threads)
    except OSError as e:
        sys.stderr.write("%s\n" % e())
        exit(2)
    sys.stderr.write("Phasing reads...\n")
    i=0
    for phasable_sample in input_data:
        sys.stderr.write('Processing alignments for sample %d (%s)...\n' % (i, phasable_sample.sample))
        header_template_file = AlignmentFile(phasable_sample.alignment_file_paths[i], 'rb')

        # Instantiate output file(s)
        output_streams = {}
        if args.output_mode == 'combined':
            output_streams["combined_output_bam"] = AlignmentFile('%s/%s.combined_phasing_results.bam' % (args.output_directory_name, phasable_sample.sample),
                                                                  'wb', template = header_template_file)
        else:
            output_streams["unphased_output_bam"] = AlignmentFile('%s/%s.unphased.bam' %  (args.output_directory_name, phasable_sample.sample),
                                                                  'wb', template = header_template_file)
            output_streams["nonphasable_output_bam"] = AlignmentFile('%s/%s.nonphasable.bam' %  (args.output_directory_name, phasable_sample.sample),
                                                                     'wb', template = header_template_file)
            if args.output_mode == 'phase_tagged':
                output_streams["phased_output_bam"] = AlignmentFile('%s/%s.phase_tagged.bam' %  (args.output_directory_name, phasable_sample.sample),
                                                                    'wb', template = header_template_file)
            elif args.output_mode == 'full':
                output_streams["maternal_output_bam"] = AlignmentFile('%s/%s.maternal.bam' %  (args.output_directory_name, phasable_sample.sample),
                                                                      'wb', template = header_template_file)
                output_streams["paternal_output_bam"] = AlignmentFile('%s/%s.paternal.bam' %  (args.output_directory_name, phasable_sample.sample),
                                                                      'wb', template = header_template_file)

        # For the negative-binomial phasing error model, we need to store phased reads for post-processing.
        phased_reads_list = []
                
        # Phase reads and write output to approriate file(s)
        for alignment in phasable_sample:
            phased_read = PhasedRead(alignment,
                                     vcf_file = phasable_sample.vcf_file_path,
                                     sample = phasable_sample.sample,
                                     powlaw_alpha = args.powlaw_alpha,
                                     powlaw_xmin = args.powlaw_xmin,
                                     multinomial_correction = args.multinomial_correction
            )

            # Make sure reads are tagged appropriately, according to output options.
            if args.add_phasing_tag_to_output:
                # Tag with estimated error rate based on fitted power-law model.
                phased_read.aligned_segment.set_tag(tag = 'ER', value = str(phased_read.phasing_error_rate), value_type='Z', replace=True)
                if phased_read.is_phased:
                    if (args.log_likelihood_threshold >= 0 and phased_read.log_likelihood_ratio >= args.log_likelihood_threshold) or phased_read.phasing_error_rate <= args.error_rate_threshold:
                        # Phased read passing the error rate/score threshold.
                        phased_read.aligned_segment.set_tag(tag = 'HP', value = str(phased_read.phase), value_type='Z', replace=True)
                    else:
                        # Phasing score does not pass threshold. Label as unphased.
                        phased_read.aligned_segment.set_tag(tag = 'HP', value = "Unphased", value_type='Z', replace=True)
                elif phased_read.is_Unphased:
                    phased_read.aligned_segment.set_tag(tag = 'HP', value = "Unphased", value_type='Z', replace=True)
                else: # Nonphasable
                    phased_read.aligned_segment.set_tag(tag = 'HP', value = "nonphasable", value_type='Z', replace=True)

            # Write output to appropriate destination file(s). All unphased/nonphasable reads
            # can be processed immediately. If not using the negative binomial error model,
            # we can process phased reads now too. For the negative-binomial error model, we
            # need to store all phased reads for later processing/relabeling/filtering.
            if args.FDR_threshold <= 0 or not phased_read.is_phased:
                write_bam_output(phased_read, output_streams, args)
            else:
                # Read is phased and we are using the negative-binomial phasing error model.
                phased_reads_list.append(phased_read)
                
        # Post-process phased reads to set FDR with the negative-binomial phasing
        # error model.
        if args.FDR_threshold > 0:
            # We need to:
            # 1) Get the average sequencing error rate from the reads (Should this be done for all reads, not just phased?)
            # 2) Get the expected number of errors (E) from the negative-binomial model
            # 3) Sort reads by LLR (Or just sort LLRs) to get lowest "acceptable" LLR (i.e., that for the E+1th value)
            # 4) Relabel reads with LLR < LLR_LIST[E+1] as unphased
            # 4b) Add tag for Z-score of some sort
            # 5) Print results
            trim_phased_reads_to_fdr_neg_binom(phased_reads_list, args.FDR_threshold, output_streams, args)
                
        # Close output file(s)(
        for outfile in output_streams.keys():
            #if outfile in locals():
            #    locals()[outfile].close()
            output_streams[outfile].close()
                
        i += 1
        
    return


def trim_phased_reads_to_fdr_neg_binom(phased_reads_list, FDR_threshold, output_streams, args):
    """
    Given a list of phased reads, control FDR at the given level
    by calculating the expected number of phasing errors based on
    a negative binomial model. Phasing decisions can be considered 
    Bernoulli trials with Psuccess = 1-(𝜀/3). Thus, the error 
    distribution may be modeled as negative binomial with N = the 
    number of phaseable reads. To control the error rate at any 
    given FDR, we use the mean of this distribution, N𝜀, as the 
    expected number of phasing errors, and reassign the N𝜀* (1-FDR)
    lowest-scoring reads as “unphased.”
    """
    r = get_average_seq_error_rate(phased_reads_list)
    p = 1 - (r / 3)
    n = len(phased_reads_list)
    E = get_expected_error_count_neg_binom(n, p)
    e = round(E * (1-FDR_threshold))
    #sys.stderr.write('p: %.5f, n: %d, E: %d, e: %d, FDR: %.5f\n' % (p, n, E, e, FDR_threshold))
    sys.stderr.write('Negative-binomial error distribution params:\nAvg. Seq. Error Rate: %.5f, Neg. Binom. p: %.5f, phased read count: %d, estimated error count: %d, FDR: %.5f\n' % (r, p, n, e, FDR_threshold))
    sorted_reads_list = sort_reads_list(phased_reads_list)
    sys.stderr.write('Reassigning %d phased reads with the lowest log-likelihood ratios as Unphased to control FDR at %.2f...\n' % (e, FDR_threshold)) 
    for i, phased_read in enumerate(sorted_reads_list):
        if i < e:
            # Relabel the first e observations as Unphased.
            phased_read.PhaseSet_max.phase = 'Unphased'
            phased_read.aligned_segment.set_tag(tag = 'HP', value = "Unphased", value_type='Z', replace=True)
        #sys.stderr.write("%s: %s, %s\n" % (phased_read.query_name, phased_read.PhaseSet_max.phase, phased_read.PhaseSet_max.max_log_likelihood_ratio))
        write_bam_output(phased_read, output_streams, args)
    

def sort_reads_list(phased_reads_list, reverse=False):
    """
    Return a sorted version of the phased_reads_list on the given attribute.
    """
    ret = sorted(phased_reads_list, key= lambda x: x.PhaseSet_max.max_log_likelihood_ratio, reverse=reverse)
    return ret


def get_average_seq_error_rate(phased_reads_list):
    """
    Get the average per-base sequencing error rate from a list of reads.
    """
    err_rates = []
    for phased_read in phased_reads_list:
        err_rates.append(phased_read.PhaseSet_max.error_rate_average_het_sites)
    return np.mean(err_rates)


def get_expected_error_count_neg_binom(n, p):
    """
    Get the expected error count from a negative-binomial distribution with
    the given parameters: n = total phased read count, p = 1-(sequencing error rate)/2.
    """
    return nbinom.stats(n, p, moments='m')


def write_bam_output(phased_read, output_streams, args):
    """
    Generic function to write a phased read as BAM output
    """
    if args.output_mode == 'combined':
        write_combined_bam_output(phased_read,
                                  output_streams["combined_output_bam"],
                                  args)
    elif args.output_mode == 'phase_tagged':
        write_phase_tagged_bam_output(phased_read,
                                      output_streams["phased_output_bam"],
                                      output_streams["unphased_output_bam"],
                                      output_streams["nonphasable_output_bam"],
                                      args)
    elif args.output_mode == 'full':
        write_full_bam_output(phased_read,
                              output_streams["maternal_output_bam"],
                              output_streams["paternal_output_bam"],
                              output_streams["unphased_output_bam"],
                              output_streams["nonphasable_output_bam"],
                              args)
    return


def write_combined_bam_output(phased_read, combined_output_bam, args):
    """
    Write phased, unphased, and nonphasable reads to a single output file.
    """
    phased_read.write_to_bam(output_bam_pysam = combined_output_bam)
    return


def write_phase_tagged_bam_output(phased_read, phased_output_bam,
                                  unphased_output_bam, nonphasable_output_bam, args):
    """
    Write phased reads to a single output file, with unphased and nonphasable
    reads written to their own files.

    """
    if phased_read.is_phased:
        if (args.log_likelihood_threshold >= 0 and phased_read.log_likelihood_ratio >= args.log_likelihood_threshold) or phased_read.phasing_error_rate <= args.error_rate_threshold:
            # Phased read passing the error-rate/score threshold.
            phased_read.write_to_bam(output_bam_pysam = phased_output_bam)
        else:
            # Phasing score does not pass threshold. Treat as unphased.
            phased_read.write_to_bam(output_bam_pysam = unphased_output_bam)
    elif phased_read.is_Unphased:
        phased_read.write_to_bam(output_bam_pysam = unphased_output_bam)
    else: # Nonphasable
        phased_read.write_to_bam(output_bam_pysam = nonphasable_output_bam)
    return
            

def write_full_bam_output(phased_read, maternal_output_bam, paternal_output_bam,
                          unphased_output_bam, nonphasable_output_bam, args):
    """
    Write maternal and paternal phased reads to separate output files, with 
    unphased and nonphasable reads written to their own files.
    """
    # TO-DO: Need to handle all possible levels of ploidy.
    if phased_read.is_phased:
        if (args.log_likelihood_threshold >= 0 and phased_read.log_likelihood_ratio >= args.log_likelihood_threshold) or phased_read.phasing_error_rate <= args.error_rate_threshold:
            if phased_read.phase == 1:
                phased_read.write_to_bam(output_bam_pysam = paternal_output_bam)
            else:
                phased_read.write_to_bam(output_bam_pysam = maternal_output_bam)
        else:
            # Phasing score does not pass threshold. Label as unphased and handle accordingly.
            phased_read.write_to_bam(output_bam_pysam = unphased_output_bam)
    elif phased_read.is_Unphased:
        phased_read.write_to_bam(output_bam_pysam = unphased_output_bam)
    else: # Nonphasable
        phased_read.write_to_bam(output_bam_pysam = nonphasable_output_bam)
    return



def phasability(args):
    """Main function for phasability mode

    Args:
        args:
    """
    
    print('phasability mode')
    #print(args)

    # First precompute error rates over a (un)reasonable range of heterozygous site counts.
    max_hets = 1000
    prior_probs = _estimate_prior_probabilities()
    # TO-DO: _error_rate_table will need to be generalized to accommodate negative-binomial FDR-setting
    _error_rate_table = _calculate_error_rate_table(max_hets, args.sequencing_error_rate,
                                                    prior_probs, args.powlaw_alpha, args.powlaw_xmin)

    # Load up the variants
    vcf_file = VariantFile(args.vcf_file_name)

    # Load up chrom_sizes/intervals
    if args.regions:
        intervals = _build_intervals_from_bed_file(args.regions)
    else:
        intervals = _build_intervals_from_chrom_sizes_dict(_chrom_sizes_dict_from_file(args.chrom_sizes))

    # Step through windows, pulling out variant counts and
    # returning the associated error rates.
    with open("%s/estimated_phasability.wig" % args.output_directory_name, 'w') as outfile:
        # Print the wiggle file header
        outfile.write('track type=wig name="Phasability" description="Phasability Plot"\n')
        for interval in intervals:
            chrom, int_start, int_end = interval

            # Make sure chromosome is in VCF
            if chrom not in list(vcf_file.header.contigs):
                continue
            
            # Print wig header
            outfile.write("fixedStep chrom=%s start=%s step=%d span=%d\n" % (chrom, int_start, args.step_size, args.step_size))

            # Step through windows within the chromosome
            win_start = int_start
            win_end = win_start + args.window_size
            while win_start < win_end and win_end <= int_end:
                n_hets = len([het for het in vcf_file.fetch(chrom, win_start, win_end)])
                if n_hets > max_hets:
                    # Have to calculate if more hets than anticipated
                    error_rate = _calculate_error_rate_from_het_count(n_hets, args.sequencing_error_rate,
                                                                      prior_probs,
                                                                      args.powlaw_alpha, args.powlaw_xmin)
                else:
                    # Use precomputed values
                    error_rate = _error_rate_table[n_hets]
                outfile.write("%.5f\n" % (1-error_rate))

                win_start += args.step_size
                win_end = win_start + args.window_size
                if win_end > int_end:
                    win_end = int_end
            
    # Clean up after ourselves
    vcf_file.close()    
    
    return


def _build_intervals_from_bed_file(bed_file):
    """
    Build a list of intervals (simply tuples of chrom, start, end) from
    a BED file.
    """
    ret = []
    with open(bed_file, 'r') as infile:
        for line in infile:
            if re.match('^#', line):  # Comment line
                continue
            rec = line.strip('\n').split()
            ret.append((rec[0], rec[1], rec[2]))
    return ret


def _build_intervals_from_chrom_sizes_dict(chrom_sizes):
    """
    Build a list of intervals from the chrom_sizes dict.
    """
    ret = []
    for chrom in chrom_sizes.keys():
        ret.append((chrom, 1, chrom_sizes[chrom]))
    return ret


def _chrom_sizes_dict_from_file(chrom_sizes):
    """
    Load chrom_sizes entries from a file into a dictionary.
    """
    chrom_sizes_dict = {}
    with open(chrom_sizes, 'r') as chrom_sizes_file:
        for line in chrom_sizes_file:
            rec = line.strip('\n').split()
            chrom_sizes_dict[rec[0]] = int(rec[1])
    return(chrom_sizes_dict)
    

def _calculate_error_rate_from_het_count(het_count, sequencing_error_rate, prior_probs, powlaw_alpha, powlaw_xmin):
    """
    Caclulate an empirical error rate given a discrete het count, sequencing
    error rate, bayesian priors, and powlaw alpha and xmin.
    """
    x = ((het_count  * np.log10(1 - sequencing_error_rate)) + max(prior_probs)) - ((het_count * np.log10(sequencing_error_rate/3)) + min(prior_probs))
    # We need a heuristic to handle the fact that values below xmin return non-probabilities.
    if x < powlaw_xmin:
        # Error rate below xmin is equivalent to a coin-toss
        return 0.5
    else:
        return powlaw_modified(x-y, powlaw_alpha, powlaw_xmin)


def _calculate_error_rate_table(max_hets, sequencing_error_rate, prior_probs, powlaw_alpha, powlaw_xmin):
    """
    Precompute empirical phasing error rates across a range of heterozygous count sites.
    Returns a list of error rates for heterozygous site counts from 0 to max_hets. These
    represent the error rates expected given a perfect match to a single phase, given a
    power law distribution with parameters powlaw_alpha and powlaw_xmin.
    """
    # This calculation is only correct for the diploid case!!

    # Numerator (all positions match)
    # For ploidy > 2, we would need to pick a specific phase and use the correct prior.
    x = np.array(range(1,max_hets+1))
    x = (x * np.log10(1 - sequencing_error_rate)) + max(prior_probs)
    
    # Denominator (all positions mismatch)
    # For ploidy > 2, this would need to be a loop over all other alleles,
    # with appropriate attention paid to using the correct prior for each.
    y = np.array(range(1,max_hets+1))
    y = (y * np.log10(sequencing_error_rate/3)) + min(prior_probs)

    # Convert to log-odds-ratios
    r = x - y

    # Extract error rates
    ret = []
    for i in range(max_hets):
        e = 0.5
        if r >= powlaw_xmin:
            e = powlaw_modified(r[i], powlaw_alpha, powlaw_xmin)
        ret.append(e)
    # Handle zero case
    ret.insert(0,0.5)
    return(ret)

def error_analysis(args):
    """Main simulation mode

    Args:
        args:
    """
    start_time = time.time()
    sys.stderr.write('\nPerforming error analysis using simulation mode...\n')
    sys.stderr.write('Analysis started %s.\n\n' % time.asctime())
    input_data = InputData(args.output_directory_path,
                           sample = args.one_sample,
                           ignore_samples = args.ignore_samples,
                           reference_sequence_input = args.reference_genome,
                           reference_sequence_input_assembly = args.reference_assembly,
                           max_threads = args.threads,
                           quiet = args.quiet_mode,
                           silent = args.silent_mode
    )
    input_data.add_haplotype_information(args.vcf_file_path)
    input_data.add_reference_sequence(args.reference_genome,
                                      reference_sequence_input_assembly = args.reference_assembly,
                                      output_directory = args.output_directory_path)

    if args.autosomal_only:
        # Only do simulations for autosomal sequences
        autosomal_ref_seq_path = create_fasta_file(
            input_fasta_file_path = args.reference_genome,
            output_fasta_file_path =  args.output_directory_path + '/reference_sequences/autosomal_' + args.reference_assembly + '_ref.fa',
            only_autosomal= True,
            chrom_sizes = args.chrom_sizes
        )
    else:
        autosomal_ref_seq_path = args.reference_genome
    
    """
    Some of what is implemented below is already contained in SimulatePhasedData.py and
    PhaseSet.py. However, I am choosing to reimplement here because I feel it will be faster
    and more efficient from a development standpoint to do so than to try and parse out how
    to use the object methods that Greg implemented. This should probably be addressed in
    a future release. --AGD
    """
    
    # Prepare the haploytpe-specific fasta
    haplotypes = [1,2]  # TO-DO: Extract the number of haplotypes from the VCF
    haplotype_specific_fasta = []
    chain_file_paths = []
    for haplotype in haplotypes:
        hsf, cfp = generate_haplotype_specific_fasta(
            haplotype,
            args.one_sample,
            autosomal_ref_seq_path,
            args.vcf_file_path,
            output_reference_sequence_path =
            args.output_directory_path +
            '/reference_sequences/autosomal_' +
            args.reference_assembly +
            '_ref_' +
            args.one_sample +
            '_hap' + str(haplotype) + '.fa',
            chain_file_path =
            args.output_directory_path +
            '/reference_sequences/autosomal_' +
            args.reference_assembly +
            '_ref_' +
            args.one_sample  +
            '_hap' + str(haplotype) + '.chain'
        )
        haplotype_specific_fasta.append(hsf)
        chain_file_paths.append(cfp)

    # Simulate reads at 1X depth for N iterations to achieve the necessary number of total reads.
    # Each 1X simulation is expected to produce ~ 190k reads.
    # TO-DO: - Simulate over all samples in VCF.
    #        - Make function args configurable through command-line options.
    hist_all = np.array([])
    hist_err = np.array([])
    bins = np.array([])
    error_lr_list = []
    error_read_list = []

    # Kludge to avoid needing a separate loop when args.simulation_min_err_obs is set:
    if args.simulation_min_err_obs > 0:
        args.simulation_depth = 9999999999999  # We just need an unreasonably large number of iterations to be safe!

    # Open a file for per-iteration simulation stats, if requested
    if args.store_simulation_data_for_iterations:
        per_iteration_stats_file = open(args.output_directory_path + "/per_iteration_simulation_error_stats.tsv", "w")
        
    sys.stderr.write("\n################ Simulating Error Data ################\n")
    for i in range(args.simulation_depth):
        sys.stderr.write("Starting iteration %d...\n" % (i+1))
        total_err_obs = len(error_lr_list)

        hist_all_i, hist_err_i, bins_i, error_lr_list_i, error_read_list_i = simulate_results_for_iteration(
            input_data,
            haplotypes,
            haplotype_specific_fasta,
            args
        )
        for j in range(len(error_lr_list_i)):
            error_lr_list.append(error_lr_list_i[j])
            error_read_list.append(error_read_list_i[j])
            
        if i > 0:
            # Combine with results of earlier iterations
            hist_all, hist_err, bins = combine_hists_and_bins(hist_all, hist_all_i,
                                                              hist_err, hist_err_i,
                                                              bins, bins_i,
                                                              args.log_likelihood_ratio_interval
            )
        else:
            hist_all, hist_err, bins = (hist_all_i, hist_err_i, bins_i)

        err_obs_i = len(error_lr_list) - total_err_obs
        sys.stderr.write("Iteration %d produced %d error observations.\n" % (i+1, err_obs_i))

        # Write cumulative error stats for each iteration to disk if requested
        if args.store_simulation_data_for_iterations:
            per_iteration_stats_file.write("iteration: %d\n" % (i+1))
            per_iteration_stats_file.write("errors for iteration: %d\tcumulative errors: %d\n" % (err_obs_i, len(error_lr_list)))
            write_simulation_stats(hist_all, hist_err, bins, per_iteration_stats_file)
            per_iteration_stats_file.write("\n")

        # Write the full error_lr_list to a file, if requested.
        if args.store_error_probs:
            # TO-DO: This file should be created once and written incrementally with
            # only reads created on this iteration.
            with open(args.output_directory_path + "/error_lr_list.txt", 'w') as outfile:
                outfile.write("%s\n" % ",".join([str(lr) for lr in error_lr_list]))

        # Write detailed information on simulated error reads, if requested.
        if args.save_read_stats:
            # TO-DO: This file should be created once and written incrementally with
            # only reads created on this iteration.
            with open(args.output_directory_path + "/error_read_stats.txt", 'w') as outfile:
                outfile.write("seq_name\tlength\tseq_error_rate_avg\tn_hets\tphase\tn_matches_to_phase\tn_mismatches_to_phase\tlog_likelihood_hap1\tlog_likelihood_hap2\tlog_likelihood_ratio\tphasing_err_rate\n")
                for err_read in error_read_list:
                    sys.stderr.write("%s\n" % (err_read.PhaseSet_max))
                    err_rate_mean = 0
                    try:
                        err_rate_mean = np.mean(err_read.PhaseSet_max.read_base_error_rates)
                    except:
                        sys.stderr.write("%s\n" % (err_read.PhaseSet_max.read_base_error_rates))
                    outfile.write("%s\t%d\t%.5f\t%d\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % (
                        err_read.read_name,
                        err_read.aligned_segment.reference_length,
                        err_rate_mean,
                        err_read.PhaseSet_max.total_hets_analyzed,
                        err_read.phase,
                        err_read.PhaseSet_max.matches[err_read.phase-1],
                        err_read.PhaseSet_max.non_matches[err_read.phase-1],
                        err_read.PhaseSet_max.log_probability_read_given_haplotype_i[0],
                        err_read.PhaseSet_max.log_probability_read_given_haplotype_i[1],
                        err_read.log_likelihood_ratio,
                        err_read.phasing_error_rate
                    ))
            
        # The rest of the kludge for using args.simulation_min_err_obs to control simulations:
        if args.simulation_min_err_obs > 0 and len(error_lr_list) >= args.simulation_min_err_obs:
            break

    # Close the per-iteration stats file if necessary.
    if args.store_simulation_data_for_iterations:
        per_iteration_stats_file.close()
            
    sys.stderr.write("\n################ Fitting Error Model ################\n")
    with open(args.output_directory_path + "/simulation_error_stats.tsv", "w") as outfile:
        write_simulation_stats(hist_all, hist_err, bins, outfile)

    # Fit the power-law model based on the observed log-likelihood-ratios and associated error rates.
    pl_fit = powerlaw.Fit(error_lr_list, method="KS")

    # Apply a linear xmin correction based on sample-size. 3.585 is the intercept term derived
    # from fitting a linear model to the difference between actual and estimated xmin values
    # over many sets of values sampled from powerlaw distributions with known parameters. We have
    # shown this correction to reduce the bias in xmin estimates across a wide range of sample sizes.
    corrected_xmin = ( pl_fit.xmin * (3.585 / np.log(len(error_lr_list))) ) + 3

    with open(args.output_directory_path + "/simulation_error_model_output.txt", 'w') as outfile:
        outfile.write("Total errors in simulations: %d\nMin log-likelihood-ratio for errors: %.3f\nMax log-likelihood ratio for errors: %.3f\nPowerlaw fitting method: %s\nNumber of observations used in powerlaw fitting: %d\nFitted powerlaw alpha: %.5f\nFitted powerlaw xmin: %.5f\nCorrected powerlaw xmin: %.5f\nFitted powerlaw minimum xmin distance (D): %.3f\nFitted powerlaw sigma: %.5f\n" % (len(error_lr_list), min(error_lr_list), max(error_lr_list), pl_fit.fit_method, pl_fit.n, pl_fit.alpha, pl_fit.xmin, corrected_xmin, pl_fit.D, pl_fit.sigma) )
        
    runtime = time.time() - start_time
    sys.stderr.write('Analysis finished %s.\nTotal runtime %.2f seconds.\n' % (time.asctime(), runtime))
    return


def write_simulation_stats(hist_all, hist_err, bins, outfile):
    # Write simulation stats to the given file.
    outfile.write("log-odds-bin\ttotal-phased\ttotal-incorrect\teffective-error-rate\n")
    for i in range(len(hist_all)):
        err_rate = -1
        if hist_all[i] > 0:
            err_rate = hist_err[i]/hist_all[i]
        outfile.write("%.5f\t%d\t%d\t%.5f\n" % (bins[i], hist_all[i], hist_err[i], err_rate))


def combine_hists_and_bins(hist_all, hist_all_i, hist_err, hist_err_i, bins, bins_i, log_likelihood_ratio_interval):
    if len(bins_i) > len(bins):
        bins = np.append(bins, np.arange(max(bins),
                                         max(bins_i)+log_likelihood_ratio_interval,
                                         log_likelihood_ratio_interval)
        )
        hist_all = sum_unequal_hists(hist_all, hist_all_i)
        hist_err = sum_unequal_hists(hist_err, hist_err_i)
    elif len(bins) > len(bins_i):
        hist_all = sum_unequal_hists(hist_all_i, hist_all)
        hist_err = sum_unequal_hists(hist_err_i, hist_err)
    else:
        hist_all = hist_all + hist_all_i
        hist_err = hist_err + hist_err_i
        
    return hist_all, hist_err, bins


def sum_unequal_hists(shorter_arr, longer_arr):
    """
    Return the sum of two numpy histograms with different maximum bins.
    """
    return longer_arr + np.append(shorter_arr, np.zeros(len(longer_arr) - len(shorter_arr)))


def simulate_results_for_iteration(input_data, haplotypes, haplotype_specific_fasta, args):
    """
    Simulate data for a single iteration and tabulate results.
    """
    hist_all = np.array([])
    hist_err = np.array([])
    bins = np.array([])
    error_lr_list = []
    error_read_list = []
    for haplotype in haplotypes:
        sys.stderr.write("Simulating reads for haplotype %d...\n" % (haplotype))
        simulated_fastq = simulate_data(haplotype_specific_fasta[haplotype-1],
                                        args.simulation_error_model,
                                        args.output_directory_path+'/simulated_reads/',
                                        args.one_sample,
                                        haplotype,
                                        quiet = args.quiet_mode,
                                        silent = args.silent_mode,
                                        no_sim = args.phase_only
        )
        try:
            sys.stderr.write("Aligning simulated reads to reference genome...\n")
            input_data.add_reads(simulated_fastq,
                                 sample = args.one_sample,
                                 reference_sequence_input = args.reference_genome,
                                 clobber=True, quiet = args.quiet_mode,
                                 silent = args.silent_mode,
                                 threads = args.threads,
                                 no_align = args.phase_only
            )
        except Exception as e:
            sys.stderr.write("%s\n" % e)
            continue
        
        # Extract log-likelihood ratios and phasing results for properly-mapped reads with >= min_het_sites
        sys.stderr.write("Phasing simulated reads...\n")
        hist_all_i, hist_err_i, bins_i, error_lr_list_i, error_read_list_i = parse_simulated_data(input_data, args)
        for j in range(len(error_lr_list_i)):
            error_lr_list.append(error_lr_list_i[j])
            error_read_list.append(error_read_list_i[j])
            
        if haplotype > 1:
            # Combine with results from other haplotype(s)
            hist_all, hist_err, bins = combine_hists_and_bins(hist_all, hist_all_i,
                                                              hist_err, hist_err_i,
                                                              bins, bins_i,
                                                              args.log_likelihood_ratio_interval
	    )
        else:
            hist_all, hist_err, bins = (hist_all_i, hist_err_i, bins_i)

        # Clean up after ourselves.
        os.remove(simulated_fastq)
        #input_data._purge_alignment_files()
        sys.stderr.write("Done.\n")
        
    return hist_all, hist_err, bins, error_lr_list, error_read_list


def simulate_data(haplotype_specific_fasta,
                  simulation_model,
                  output_directory,
                  sample,
                  haplotype,
                  depth = 1,
                  difference_ratio = '23:31:46',
                  length_mean = 25000,
                  length_max = 1000000,
                  length_min = 100,
                  length_sd = 20000,
                  accuracy_min = 0.01,
                  accuracy_max = 1.00,
                  accuracy_mean = 0.80,
                  quiet = False, silent = False,
                  no_sim = False):
    """
    Simulate reads for a single sample and haplotype.
    TO-DO: Make all optional args accessible as command-line options.
    """
    simulated_fastq = simulate_reads_pbsim2(
        reference_sequence = haplotype_specific_fasta,
        depth = depth,
        simulation_mode = simulation_model,
        difference_ratio = difference_ratio,
        length_mean = length_mean,
        length_max = length_max,
        length_min = length_min,
        length_sd = length_sd,
        accuracy_min = accuracy_min,
        accuracy_max = accuracy_max,
        accuracy_mean = accuracy_mean,
        output_directory = output_directory,
        sample = sample,
        haplotype = haplotype,
        quiet = quiet,
        silent = silent,
        no_sim = no_sim
    )
    
    return simulated_fastq


def parse_simulated_data(input_data, args):
    """
    Takes simulated alighments and parses out log-likelihood ratios and associated error counts.
    """
    log_ratios = []
    is_phased_correctly = []
    error_lr_list = []
    error_read_list = []
    for phasable_sample in input_data:
        # To-Do: Work in ability to select/ignore specific samples
        for alignment in phasable_sample:
            phased_read = PhasedRead(alignment, vcf_file = phasable_sample.vcf_file_path, sample = 'HG001', evaluate_true_alignment = True)
            if phased_read.alignment_is_mapped and phased_read.matches_true_alignment and phased_read.is_phased and int(phased_read.aligned_segment.get_tag('HS')) > 1:
                log_ratios.append(phased_read.log_likelihood_ratio)
                is_phased_correctly.append(phased_read.is_phased_correctly)
                if not phased_read.is_phased_correctly:
                    #sys.stderr.write("%s\n%s\n" % (phased_read, phased_read.is_phased_correctly))
                    error_lr_list.append(phased_read.log_likelihood_ratio)
                    error_read_list.append(phased_read)
        
    # Tabulate the total counts and number of errors at each log-odds level (binned in steps of bin_width).
    log_ratios = np.array(log_ratios)
    is_phased_correctly = np.array(is_phased_correctly, dtype=bool)
    bins = np.arange(0, np.ceil(max(log_ratios))+args.log_likelihood_ratio_interval, args.log_likelihood_ratio_interval)
    hist_all, edges = np.histogram(log_ratios, bins)
    hist_err, edges = np.histogram(log_ratios[np.nonzero(np.where(is_phased_correctly, is_phased_correctly == False, 1))], bins)
    
    return hist_all, hist_err, bins, error_lr_list, error_read_list


def getArgs() -> object:
    """Parses arguments:"""
    ################################################################################################################
    ############## LRphase Arguments ###############################################################################
    
    lrphase_parser: argparse.ArgumentParser = argparse.ArgumentParser(
        prog = 'LRphase', description = 'Tools for phasing individual long reads using haplotype information.'
    )
    
    lrphase_parser.add_argument('--version', action = 'version', version = __version__)

    lrphase_subparsers: argparse._SubParsersAction = lrphase_parser.add_subparsers(
        title = '[LRphase modes]', dest = 'mode', description = 'Choose which mode to run:',
        help = 'mode must be added as first argument (ex: LRphase phasing)', required = True
    )

    
    ################################################################################################################
    ############## Phasing Mode Arguments ##########################################################################
    
    phasing_parser: argparse.ArgumentParser = lrphase_subparsers.add_parser(
        'phasing', description = "Tool for phasing individual long reads using haplotype information."
    )
    
    ############## Phasing Mode Required Arguments ##############
    phasing_parser_required: argparse._ArgumentGroup = phasing_parser.add_argument_group('Required', 'Required for phasing')
    
    phasing_parser_required.add_argument(
        '-v', '--vcf', required = True,
        help = 'Path to vcf file with haplotype information that will be used for phasing. (Must be in .vcf.gz format '
        'with tabix index in same folder. If .vcf file is provided, bgzip and tabix must be installed and '
        'available on PATH because LRphase will attempt to convert it.  EX: -v GM12878_haplotype.vcf.gz)',
        dest = 'vcf_file_name', metavar = '<VCF_FILE>'
    )
    
    phasing_parser_required.add_argument(
        '-i', '--input_reads', required = True,
        help = 'Path to sequencing file (.fasta) or alignment file (.bam or .sam) of long reads that will be used for phasing. If either a .sam file is provided or an index is not found, .sam and .bam file will be sorted and indexed with SAMtools. Sorted.bam files should be in same directory as their index (.sorted.bam.bai). EX: -a data/minion_GM12878_run3.sorted.bam, -i minion_GM12878_run3.sam) Path to long read file in .fastq format that will be used for alignment and phasing (ex: -i minion_GM12878_run3.fastq). **** NOTE: the -r/--reference argument is REQUIRED if using input in fastq format! ****',
        #dest = 'long_read_inputs', metavar = 'long-read input file', action = 'append', nargs = '*'
        dest = 'long_read_inputs', metavar = '<SAM/BAM/FASTQ>'
    )

    ############## Phasing Mode Optional Arguments ##############
    phasing_parser_optional: argparse._ArgumentGroup = phasing_parser.add_argument_group('Optional', 'Useful, but (mostly) not required for phasing.')

    phasing_parser_optional.add_argument(
        '-o', '--output_directory_name', required = False,
        default = '.',
        help = 'Name given to directory where results will be output (ex: -o minion_GM12878_run3_phasing_output)',
        dest = 'output_directory_name', metavar = "</path/to/output>"
    )
    
    phasing_parser_optional.add_argument(
        '-r', '--reference', required = False,
        help = 'Path to reference genome sequence file. REQUIRED if -i is used to specify reads in fastq format to be '
        'aligned prior to phasing. (file types allowed: .fa, .fna, fasta. EX: -r data/reference_hg38.fna)',
        dest = 'reference_genome', metavar = '<REF_FASTA>'
    )
    
    phasing_parser_optional.add_argument(
        '-A', '--reference_assembly', required = False, default="hg38",
        help = 'Assembly for the reference genome. EX: -A hg38.',
        dest = 'reference_assembly', metavar = '<ASSEMBLY_NAME>'
    )

    phasing_parser_optional.add_argument(
	'-t', '--threads', type=int, required=False, default=3,
        help = 'Number of threads to use for mapping and indexing steps.',
        dest = 'threads', metavar = '<THREADS>'
    )
    phasing_parser_optional.add_argument(
        '-q', '--quiet', help = 'Output to stderr from subprocesses will be muted.', action = 'store_true', dest = 'quiet_mode'
    )

    phasing_parser_optional.add_argument(
        '-S', '--silent', help = 'Output to stderr and stdout from subprocesses will be muted.', action = 'store_true', dest = 'silent_mode'
    )
    
    ############## Phasing Mode Output Options ##############
    phasing_parser_output: argparse._ArgumentGroup = phasing_parser.add_argument_group(
        'Output options', 'Options for writing output to BAM file(s).'
    )

    phasing_parser_output.add_argument(
        '-H', '--omit_phasing_tag',
        #help = 'Do not tag reads with "HP:i:1" (maternal) and "HP:i:2" (paternal) tags to indicate phase assignments.',
        help=argparse.SUPPRESS,
        dest = 'add_phasing_tag_to_output', action = 'store_false'
    )
    
    phasing_parser_output.add_argument(
        '-P', '--omit_phase_set_tag',
        #help = 'Do not tag reads with "PS:i:x" tags, which label reads according to the phase set that was indicated in the vcf record used to assign a read to a phase.',
        help=argparse.SUPPRESS,
        dest = 'add_phase_set_tag_to_ouput', action = 'store_false'
    )

    phasing_parser_output.add_argument(
        '-Q', '--omit_phasing_quality_tag',
        #help = 'Do not tag reads with "PC:i:x" tags, where x is an estimate of the accuracy of the phasing assignment in phred-scale.',
        help=argparse.SUPPRESS,
        dest = 'add_phasing_quality_tag_to_output', action = 'store_false'
    )

    phasing_parser_output.add_argument(
        '-O', '--output_mode', type=str, required=False, default="combined",
        choices=['combined', 'phase_tagged', 'full'],
        help = 'Specify whether/how phased, unphased, and nonphasable reads are printed to output. Modes available:\n\tcombined: All reads will be written to a common output file. The phasing tag (HP:i:N) can be used to extract maternal/paternal phased reads, unphased reads, and nonphasable reads.\n\tphase_tagged: Phased reads for both maternal and paternal phases will be written to a single output file, while unphased and nonphasable reads will be written to their own respective output files.\n\tfull: Maternal, paternal, unphased, and nonphasable reads will be printed to separate output files.',
        action = 'store'
    )

    ############## Multiple sample handling/phase set options for phasing mode ##############
    #phasing_parser_multiple_sample = phasing_parser.add_argument_group(
    #    'Sample options',
    #    'Phasing options for files containing multiple samples or haplotypes. By default the input files are assumed to belong to the same haplotype and from the same sample. '
    #)
    
    #phasing_parser_multiple_sample.add_argument(
    phasing_parser_output.add_argument(
        '-s', '--one_sample', required = False,
        #help = 'Use the --one_sample option to phase a specific sample present in the input reads and vcf file. (-s HG001)',
        help=argparse.SUPPRESS,
        metavar = '<SAMPLE_NAME>'
    )
    
    #phasing_parser_multiple_sample.add_argument(
    phasing_parser_output.add_argument(
        '--ignore_samples',
        #help = 'Use the --ignore_samples option to ignore sample labels. The first sample column in the VCF will be used and reads will not be matched using RG tags, samples, or phase sets.',
        help=argparse.SUPPRESS,
        action = 'store_true'
    )
    
    ############## Statistical and error options for phasing mode ##############
    phasing_parser_stats_error: argparse._ArgumentGroup = phasing_parser.add_argument_group(
        'Statistical options for phasing model',
        'Options to modify thresholds and error parameters involved in phasing decisions.'
    )

    phasing_parser_stats_error.add_argument(
        '-F', '--FDR_threshold',
        type = float,
        required = False,
        default = 0.05,
        help = 'Control the false discovery rate at the given value using a negative-binomial estimate of the number of phasing errors (N) given the average per-base sequencing error rate observed among all phaseable reads. Phased reads are sorted by their observed log-likelihood ratios and the bottom N*(1-FDR) reads will be reassigned to the "Unphased" set. Set this to zero to skip this step and return all phasing predictions.'
    )

    phasing_parser_stats_error.add_argument(
        # Until the generalized multinomial coefficient is implemented, error models
        # other than the default (0) will not produce valid probabilities. 
        '--sequencing_error_model',
        required = False,
        default = 0,
        choices = [0,1,2],
        #help = 'Use this option to choose how to estimate sequencing error rates: 0: estimate per-base error rate as an average per read. 1: estimate per-base error rate locally around each het site. 2: Calculate per-base error rate using base quality scores (WARNING: do not use option 2 unless you are sure that the basecaller reported actual error rates).',
        help=argparse.SUPPRESS,
        dest = 'error_model',
        type = int
    )

    # Options related to the powerlaw error rate approximation are hidden with
    # defaults set to disable filtering of results based on powerlaw error rate
    # estimates. -E 1.0 stipulates that no predictions will be tossed out based on
    # empirical error rates, which should not exceed 0.5.
    phasing_parser_stats_error.add_argument(
        '-E', '--error_rate_threshold',
        type = float,
        required = False,
        #default = 0.01,
        default = 1.0,
        #help = 'Error rate threshold on phasing results. This threshold equates to the estimated phasing error rate for an experiment, such that a threshold of 0.05 should be equivalent to a 5%% false-discovery rate. This rate is estimated based on a fitted power-law distribution with alpha and xmin parameters supplied with the --powlaw_alpha and --powlaw_xmin options. These parameters may be estimated by running the LRphase error_analysis mode.',
        help=argparse.SUPPRESS,
        dest = 'error_rate_threshold', metavar = '<MAX_ERROR_RATE>'
    )

    phasing_parser_stats_error.add_argument(
        '--log_likelihood_threshold',
        type = float,
        required = False,
        default = -1.0,
        help = 'Use a hard threshold on log-likelihood ratios when phasing reads. Results will only be printed for predicted phasings with log-likelihood ratios equal to or greater than this threshold. Setting this to zero will cause all reads to be assigned to the phase to which they share the greatest number matches. Log-likelihood ratios will still be reported in the output in this case, but are not used for phasing decisions.',
        metavar = '<MIN_LIKELIHOOD_RATIO>'
    )

    phasing_parser_stats_error.add_argument(
        #This option is permanently disabled!
        #When using sequencing error models other than the default (shared error rate 
        #across all positions in a read), we need this coefficient to produce valid
        #probabilities. When local error rates are being used for each het site, we
        #still need a coefficient, but must calculate it differently, so this option
        #should never be used!
        # Note that we are storing False here!
        '--no_multcoeff', required = False, default = False, action = 'store_true',
        #help = 'Do not apply the multinomial coefficient in the likelihood calculation. Default=False (The multinomal coefficient will be used.)',
        help=argparse.SUPPRESS,
        dest = 'multinomial_correction'
    )

    phasing_parser_stats_error.add_argument(
        '--powlaw_alpha', type = float, required = False, default = 5.4,
        #help = 'Alpha parameter for the fitted power law distribition used in error-rate calculations. This can be obtained by running LRphase error_analysis mode. Default = 5.4',
        help=argparse.SUPPRESS,
        dest = 'powlaw_alpha', metavar="<ALPHA>"
    )

    phasing_parser_stats_error.add_argument(
        '--powlaw_xmin', type = float, required = False, default = 4.3,
        #help = 'Xmin parameter for the fitted power law distrubution used in error-rate calculations. This can be obtained by running LRphase error_analysis mode. Default = 4.3',
        help=argparse.SUPPRESS,
        dest = 'powlaw_xmin', metavar = '<XMIN>'
    )

    phasing_parser.set_defaults(func = phasing)

    
    """
    # Since it is not clear how to evaluate the empirical phasing error rate without the powerlaw
    # approximation, which has proven difficult to reliably attain, this mode is not discussed in
    # the 2022 manuscript. Therefore, this mode is disabled in the current (1.X) released versions.
    #########################################################################################################################
    ############## Phasability Mode Arguments ###############################################################################
    
    phasability_parser = lrphase_subparsers.add_parser(
        'phasability',
        description = 'Tool that uses haplotype information as input and outputs predictions for how well LRphase will perform for a sequencing dataset with a given N50 and sequencing error rate. Phasability is defined as the probability of correctly phasing a read of width W that matches all variants in a single phase, calculated as (1 - empirical_phasing_error_rate) for any window of width W. Overlapping windows of width W are slid along the genome at intervals of step size S and results are reported in fixed-step wig format. Individual bins in the wig file describe phasability for the W-width window starting at the leftmost position of the bin.  Can be used to evaluate phasability genome-wide or at optional regions.'
    )
    
    ############## Phasability Mode Required Arguments ##############
    phasability_parser_required: argparse._ArgumentGroup = phasability_parser.add_argument_group(
        'Required', 'Required for phasability analysis'
    )
    
    phasability_parser_required.add_argument(
        '-o', '--output_directory_name', required = True,
        help = 'Name given to directory where results will be output (ex: -o GM12878_phasability)',
        dest = 'output_directory_name', metavar = '</path/to/output>'
    )
    
    phasability_parser_required.add_argument(
        '-v', '--vcf', required = True,
        help = 'Path to vcf file with haplotype information that will be used for phasability analysis. (Must be in '
        '.vcf.gz format with tabix index in same folder. If .vcf file is provided, bgzip and tabix must be '
        'installed and available on PATH because LRphase will attempt to convert it. EX: -v '
        'GM12878_haplotype.vcf.gz)',
        dest = 'vcf_file_name', metavar = '<VCF_FILE>'
    )

    phasability_parser_required.add_argument(
        '-C', '--chrom_sizes', required = True, type=str,
        help = 'Chromosome sizes file. Required when using --autosomal_only',
        dest = 'chrom_sizes', metavar = '<CHROM_SIZES_FILE>'
    )
    
    ############## Phasability Mode Optional Arguments ##############
    phasability_parser_optional = phasability_parser.add_argument_group('Optional', 'Optional for phasability mode')
    
    phasability_parser_optional.add_argument(
        '-w', '--window_size', required = False, type=int, default = 25000,
        help = 'Window size for calculating phasability. Set this to a value near the N50 of the sequencing experiment you are evaluating. (default = 25000bp)',
        dest = 'window_size', metavar = '<N>'
    )
    
    phasability_parser_optional.add_argument(
        '-t', '--step_size', required = False, type = int, default = 1000,
        help = 'Distance between start positions of the overlapping windows used for phasability analysis (default = 1000bp)',
        dest = 'step_size', metavar = '<N>'
    )

    phasability_parser_optional.add_argument(
        '-e', '--sequencing_error_rate', required = False, type=float, default = 0.1,
        help = 'Estimaated per-base sequencing error rate. Typically ~0.1 for Nanopore sequencing. (Default = 0.1)',
        dest = 'sequencing_error_rate', metavar = '<ERROR_RATE>'
    )
    
    phasability_parser_optional.add_argument(
        '--powlaw_alpha', type = float, required = False, default = 4.5,
        help = 'Alpha parameter for the fitted power law distribition used in error-rate calculations. This can be obtained by running LRphase error_analysis mode. Default = 4.5',
        dest = 'powlaw_alpha', metavar="<ALPHA>"
    )

    phasability_parser_optional.add_argument(
        '--powlaw_xmin', type = float, required = False, default = 2.0,
        help = 'Xmin parameter for the fitted power law distrubution used in error-rate calculations. This can be obtained by running LRphase error_analysis mode. Default = 2.0',
        dest = 'powlaw_xmin', metavar = '<XMIN>'
    )
    
    phasability_parser_optional.add_argument(
        '-R', '--regions', required = False, type=str,
        help = 'Regions in which to estimate phasability, in BED format.',
        dest = 'regions', metavar = '<BED_FILE>'
    )
        
    phasability_parser_optional.add_argument(
        '-q', '--quiet', help = 'Output to stderr from subprocesses will be muted.', action = 'store_true', dest = 'quiet_mode'
    )

    phasability_parser_optional.add_argument(
	'-S', '--silent', help = 'Output to stderr and stdout from subprocesses will be muted.', action = 'store_true', dest = 'silent_mode'
    )
    
    phasability_parser.set_defaults(func = phasability)
    """

    """
    # This mode is disabled in the 1.X releases, published with the 2022 manuscript.
    ##################################################################################################################
    ############## error_analysis Mode Arguments #####################################################################
    
    error_analysis_parser = lrphase_subparsers.add_parser(
        'error_analysis',
        description = 'Tool for estimating error rates in a dataset given a (set of) haplotype(s) and a reference genome. Simulated reads with known phase are generated based on the inputs, and processed through the LRphase phasing mode. Results are used to estimate parameters for the powerlaw distribution modeling empirical error rates within the data. (See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085777, and http://arxiv.org/abs/0706.1062 for more information.). Results of powerlaw model fittign are written to "simulation_error_model_output.txt" and a text report of error rates across the observed range of log-likelihood ratios is written to "simulation_error_stats.tsv" in the output directory. These files can be used to obtain the optimal parameters for setting thresholds for phasing mode: To use the fitted powerlaw model to set an error-rate based threshold (recommended), the "Fitted powerlaw alpha" and "Corrected powerlaw xmin" within "simulation_error_model_output.txt" represent the recommended values for the --powlaw_alpha and --powlaw_xmin options, respectively. This file also reports "Fitted powerlaw xmin", which is the actual xmin estimate from model fitting. This estimate tends to be systematically low, especially with fewer than 25,000 error observations, thus "Corrected powerlaw xmin" has been adjusted using a linear correction based on sample-size. If you are simulating more than 25,000 observations, you may wish to use the fitted value rather than the corrected value. To set a hard threshold on log-odds scores (NOT recommended), the values within "simulation_error_stats.tsv" can be used to build an eCDF from which a log-odds threshold can be chosen at any desired empirical error rate, which can be supplied to phasing mode using the --error_rate_threshold option.'
    )
    
    ############## error_analysis Mode Required Arguments ##############
    
    error_analysis_parser_required: argparse._ArgumentGroup = error_analysis_parser.add_argument_group(
        'Required', 'Required for error_analysis'
    )
    
    error_analysis_parser_required.add_argument(
        '-o', '--output_directory_path', required = True,
        help = 'Name given to directory where results will be output (ex: -o minion_GM12878_run3_phasing_output)',
        dest = 'output_directory_path', metavar = '</path/to/output_directory>'
    )
    
    error_analysis_parser_required.add_argument(
        '-v', '--vcf', required = True,
        help = 'Path to vcf file with haplotype information that will be used for error_analysis. (Must be in .vcf.gz format with tabix index in same folder. If .vcf file is provided, bgzip and tabix must be installed and available on PATH because LRphase will attempt to convert it.  EX: -v GM12878_haplotype.vcf.gz)',
        dest = 'vcf_file_path', metavar = '<VCF_FILE>'
    )

    error_analysis_parser_required.add_argument(
        '-p', '--simulation_error_model', required = True,
        help = 'Path to the error model for simulated sequence generation',
        dest = 'simulation_error_model', metavar = '<ERROR_MODEL_FILE>'
    )
    
    error_analysis_parser_required.add_argument(
        '-r', '--reference', required = True,
        help = 'Path to reference genome sequence file.',
        dest = 'reference_genome', metavar = '<GENOME_FASTA>'
    )

    error_analysis_parser_required.add_argument(
        '-A', '--reference_assembly', required = False, default='hg38', type=str,
        help = 'Reference genome assembly, e.g., hg38.',
        dest = 'reference_assembly', metavar = '<ASSEMBLY_NAME>'
    )

    error_analysis_parser_required.add_argument(
        '-t', '--threads', type=int, required=False, default=3,
        help = 'Number of threads to use for mapping step.',
        dest = 'threads', metavar='<N>'
    )
    
    ############## Additional methods for evaluating phasing accuracy using simulated reads  ##############
    error_analysis_parser_simulation: argparse._ArgumentGroup = error_analysis_parser.add_argument_group(
        'Evaluation of phasing using simulated reads',
        'Simulated reads are used to estimate the accuracy of phasing decisions and assign quality scores.'
    )

    error_analysis_parser_simulation.add_argument(
        '--ignore_phase_set',
        help = 'Compare the quality scores in a phased bam to error rates determined from simulated reads. All input '
        'must be phased bam files that  -i for ath to phased with phased reads to output results from '
        'simulated analyses.',
        action = 'store_true'
    )

    error_analysis_parser_simulation.add_argument(
        '-l', '--log_likelihood_ratio_interval', required = False, type=float, default = 0.5,
        help = 'Bin width for tracking phasing error at each log-likelihood-ratio level.',
        dest = 'log_likelihood_ratio_interval', metavar = '<N>'

    )

    error_analysis_parser_simulation.add_argument(
	'--simulation_depth', required = False, default=25, type=int,
        help = 'Total read depth coverage for all simulations.',
        dest = 'simulation_depth', metavar = '<N>'
    )

    error_analysis_parser_simulation.add_argument(
	'--simulation_min_err_obs', required = False, default=-1, type=int,
        help = 'Minimum number of error observations to require before stopping simulations. It is recommended that this be set >= 1000 for reliable estimation of the xmin parameter, though as few as 50 observations is sufficient to reliably estimate alpha if xmin is already known.',
        dest = 'simulation_min_err_obs', metavar = '<N>'
    )

    error_analysis_parser_simulation.add_argument(
        '--sample', required = False, type=str,
        help = 'Sample name in VCF to use for haploytpe simulation.',
        dest = 'one_sample', metavar = '<SAMPLE_NAME>'
    )

    error_analysis_parser_simulation.add_argument(
        '--ignore_samples', required = False,
        help = 'Ignore VCF samples. Simulations will be based on all samples.',
        dest = 'ignore_samples', action='store_true'
    )

    error_analysis_parser_simulation.add_argument(
        '--autosomal_only', required = False,
        help = 'Only simulate sequences on autosomes.',
        dest = 'autosomal_only', action='store_true'
    )

    error_analysis_parser_simulation.add_argument(
        '--chrom_sizes', required = False, type=str,
        help = 'Chromosome sizes file. Required when using --autosomal_only',
        dest = 'chrom_sizes', metavar = '<FILE>'
    )

    error_analysis_parser_simulation.add_argument(
        '--store_error_probs', required = False,
        help = 'Write the list of likelihood ratios for all error observations to a comma-delimited text file.',
        dest = 'store_error_probs', action='store_true'
    )

    error_analysis_parser_required.add_argument(
        '--save_detailed_read_stats', default=False,
        help = 'Save detailed information on simulated error reads. Data will be stored to <output_dir>/error_read_stats.txt',
        dest = 'save_read_stats', action = 'store_true'
    )

    error_analysis_parser_required.add_argument(
        '--store_simulation_data_for_iterations', default=False,
        help = 'Save simulated error stats for each iteration to disk. Data will be stored to <output_dir>/per_iteration_simulation_error_stats.tsv',
        dest = 'store_simulation_data_for_iterations', action = 'store_true'
    )

    error_analysis_parser_simulation.add_argument(
	'-q', '--quiet', help = 'Output to stderr from subprocesses will be muted.', action = 'store_true', dest = 'quiet_mode'
    )

    error_analysis_parser_simulation.add_argument(
	'-S', '--silent',
        help = 'Output to stderr and stdout from subprocesses will be muted.',
        action = 'store_true', dest = 'silent_mode'
    )
    
    # For debugging output
    error_analysis_parser_simulation.add_argument(
        '--phase_only',
        #help = 'Use existing simulated alignment data to perform phasing step only.',
        help=argparse.SUPPRESS,
        action = 'store_true'
    )
    
    error_analysis_parser.set_defaults(func = error_analysis)
    """
    
    #####################################################################################################################
    ############## Parse arguments ######################################################################################    
    args = lrphase_parser.parse_args()
    args.func(args)
    return args

#@profile
def main():
    # This only needs to call the getArgs method, which will then dispatch functions for the indicated runtime mode.
    args: argparse.Namespace = getArgs()
    
    
if __name__ == '__main__':
    main()

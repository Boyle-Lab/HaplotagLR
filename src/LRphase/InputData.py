# coding=utf-8
import sys, os
import random
import subprocess
import time
from collections import defaultdict
from datetime import date
from shutil import copyfileobj
from typing import Any, Dict, Tuple, List, Optional, Union

import pysam
import requests
from LRphase import PhasableSample
from LRphase.PhasedRead import *
from pysam import VariantFile
from LRphase import urls


def _prepare_output_directory(output_directory_path: str) -> object:#Optional[str]:
    """

    Args:
        output_directory_path (object): 
    """
    if os.path.isfile(output_directory_path):
        sys.stderr.write('Output directory exists as a file. Use -o to specify the name to be given to the output folder. May also provide a relative path with a name to create output directory in a different location (EX: -o path/to/name). Do not specify a file.\n')
        return
    elif os.path.exists(output_directory_path):
        sys.stderr.write('%s already exists, WARNING: New results will overwrite old files.\n' % os.path.abspath(output_directory_path))
    elif not os.path.exists(output_directory_path):
        os.mkdir(output_directory_path)
        sys.stderr.write("%s created\n" % os.path.abspath(output_directory_path))
    output_directory = os.path.abspath(output_directory_path)
    return output_directory


def _prepare_summary_file(output_directory):
    sys.stderr.write('\n############## Preparing summary_file ##############\n')
    
    summary_file_path = os.path.abspath('%s/%s' % (output_directory, 'summary.txt'))
    if os.path.exists(summary_file_path):
        sys.stderr.write('%s already exists, WARNING: New results will overwrite old files.\n' % os.path.abspath(summary_file_path))
    elif not os.path.exists(summary_file_path):
        sys.stderr.write('%s created.\n' % os.path.abspath(summary_file_path))
    summary_file = open(summary_file_path, 'w')
    summary_file.write('Start date: %s\n' % str(date.today()))
    start_time = time.time()
    sys.stderr.write('%s created.\n' % summary_file_path)
    
    return summary_file_path


def _pair_sample_with_vcf(sample, sample_to_vcf_file_dict, ignore_phase_sets):
    if str(sample) in sample_to_vcf_file_dict.keys():
        for vcf_file_path in sample_to_vcf_file_dict[str(sample)].items():
            if ignore_phase_sets:
                return vcf_file_path[0], True
            else:
                return vcf_file_path[0], vcf_file_path[1]
    else:
        return


def _sample_to_alignment_files(sample_to_vcf_file_dict: dict, RG_ID_dict: dict) -> object:
    """
    Create a dictionary relating sample names to alignment files.
    """
    sample_to_alignment_files = {}
    for sample in sample_to_vcf_file_dict:
        #sys.stderr.write("%s\n" % sample)
        sample_to_alignment_files[sample] = {}
        for file in RG_ID_dict.keys():
            #sys.stderr.write("%s\n" % file)
            if len([file for key in RG_ID_dict[file].keys() if RG_ID_dict[file][key]['SM'] == sample]) > 0:
                #sys.stderr.write("%s\n" % file)
                for pair in [{
                        alignment_file_path:any(
                                [RG_ID_dict[alignment_file_path][ID]['RG_tags'] for ID in
                                 RG_ID_dict[alignment_file_path].keys()]
                                )
                        } for alignment_file_path in RG_ID_dict.keys() if any(
                        [RG_ID_dict[alignment_file_path][ID]['SM'] == sample for ID in
                         RG_ID_dict[alignment_file_path].keys()]
                        )]:
                    #sys.stderr.write("%s\n" % pair)
                    sample_to_alignment_files[sample][list(pair.keys())[0]] = list(pair.values())[0]
    return sample_to_alignment_files


def _align_long_reads_fastq(long_reads_fastq_path, reference_sequence_input, output_directory, threads, quiet = False, silent = False, no_align = False):
    """
    Align specified reads file to reference genome via minimap2.
    """
    if not reference_sequence_input:
        sys.stderr.write('no reference genome was provided for %s. This file will be skipped.\n' % long_reads_fastq_path)
        return

    if not no_align:
        sys.stderr.write('\n################ Beginning Sequence Alignment ################\n')
    
    start_process_time = time.time()
    reference_sequence_input_path = os.path.abspath(reference_sequence_input)
    if str(long_reads_fastq_path).lower().endswith('.fastq'):
        sam_alignment_path = os.path.abspath(
                '%s/%s%s' % (
                        output_directory, os.path.splitext(os.path.basename(long_reads_fastq_path))[0],
                        '_alignment.sam')
                )
    elif str(long_reads_fastq_path).lower().endswith('fastq.gz'):
        sam_alignment_path = os.path.abspath(
                '%s/%s%s' % (
                        output_directory,
                        os.path.splitext(os.path.splitext(os.path.basename(long_reads_fastq_path))[0])[0],
                        '_alignment.sam')
                )

    mapping_cmd_args = ['minimap2', '-ax', 'map-ont', '-Y', '-L', '--secondary=no',
                        '--MD', reference_sequence_input_path,
                        long_reads_fastq_path,
                        '-o', sam_alignment_path,
                        '-t', str(threads)]
    if not no_align:
        if silent:
            subprocess.run(mapping_cmd_args, check = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
        elif quiet:
            subprocess.run(mapping_cmd_args, check = True, stderr = subprocess.DEVNULL)
        else:
            subprocess.run(mapping_cmd_args, check = True)
    
    end_process_time = time.time()
    total_process_time = end_process_time - start_process_time
    
    sys.stderr.write('Alignment finished in %.2f seconds.\n' % total_process_time)
    
    return sam_alignment_path


def _sort_and_index_alignment_file(long_reads_alignment_path, output_directory: object, threads: int = 1):
    start_process_time = time.time()
    
    sorted_bam_file_path = os.path.abspath(
            '%s/%s%s' % (
                    output_directory, os.path.splitext(os.path.basename(long_reads_alignment_path))[0],
                    '_sorted.bam')
            )
    sys.stderr.write('\n############## Sorting and indexing bam file ##############\n')
    #pysam.sort('-O', 'BAM', long_reads_alignment_path, '-o', sorted_bam_file_path, '-@', str(threads))
    # TO-DO: Add siient/quiet args for subprocess command output
    sort_cmd_args = ['samtools', 'sort',
                    '-O', 'BAM',
                    '-o', sorted_bam_file_path,
                    '-@', str(threads),
                    long_reads_alignment_path]
    subprocess.run(sort_cmd_args, check = True)
    sys.stderr.write('Created %s.\n' % sorted_bam_file_path)
    #pysam.index('-@', str(threads), sorted_bam_file_path)
    index_cmd_args = ['samtools', 'index',
                      '-@', str(threads),
                      sorted_bam_file_path]
    subprocess.run(index_cmd_args, check = True)
    sys.stderr.write('Created %s.bai.\n' % sorted_bam_file_path)
    os.remove(long_reads_alignment_path)
    end_process_time = time.time()
    total_process_time = end_process_time - start_process_time
    sys.stderr.write('Sorting and indexing finished in %.2f seconds.\n' % total_process_time)
    return sorted_bam_file_path


def _prepare_alignment(output_directory: str, long_reads_alignment_path: str, threads: int = 1) -> str:
    start_process_time = time.time()

    # Check integrity of alignment file
    cmd = "samtools quickcheck " + long_reads_alignment_path
    exit_status = os.system(cmd)
    if exit_status != 0:
        # Bam file could not be read.
        raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_prepare_alignment_177)\n" % (long_reads_alignment_path))
        #pass

    long_read_file_pysam = pysam.AlignmentFile(long_reads_alignment_path)
    if long_read_file_pysam.format == 'BAM':
        if long_read_file_pysam.has_index():
            sorted_bam_file_path = os.path.abspath(long_reads_alignment_path)
            sys.stderr.write('%s is a valid alignment file with an index.\n' % long_reads_alignment_path)
        elif not long_read_file_pysam.has_index():
            sys.stderr.write('%s is a .bam file but the index cannot be found. Sorting and indexing bam file.\n' % long_reads_alignment_path)
            sorted_bam_file_path = _sort_and_index_alignment_file(
                long_reads_alignment_path,
                output_directory,
                threads=threads
            )
    elif long_read_file_pysam.format == 'SAM':
        sys.stderr.write('%s is a .sam file. Converting to binary (.bam), sorting, and indexing bam file.\n' % long_reads_alignment_path)
        sorted_bam_file_path = _sort_and_index_alignment_file(
            long_reads_alignment_path,
            output_directory,
            threads=threads
        )
    else:
        sys.stderr.write('Error: Pysam does not recognize %s as being in SAM or BAM format. If aligned reads are provided as input they must be in proper .sam or .bam format.' % long_reads_alignment_path)
        raise Exception('Improper alignment file format.')
    return sorted_bam_file_path


def _unique_RG_IDs_from_RG_tags(RG_ID_dict: dict, unique_RG_IDs: dict, alignment_file_path: str) -> object:
    """

    Args:
        RG_ID_dict (dict):
        unique_RG_IDs (dict):
        alignment_file_path (str):
    """
    # Check integrity of alignment file
    cmd = "samtools quickcheck " + long_reads_alignment_path
    exit_status = os.system(cmd)
    if exit_status != 0:
        # Bam file could not be read.
        raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_unique_RG_217)\n" % (long_reads_alignment_path))
        #pass
    
    with pysam.AlignmentFile(alignment_file_path, 'rb') as bam_file:
        RG_tags = bam_file.header.get('RG')
    if RG_tags is None:
        RG_ID_dict[str(alignment_file_path)] = 'No RG tags'
    else:
        RG_ID_dict[str(alignment_file_path)] = {}
        for RG_tag in RG_tags:
            RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])] = {}
            RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['DS'] = str(RG_tag['DS'])
            RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['SM'] = str(RG_tag['SM'])
            RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['RG_tags'] = True
            if str(RG_tag['ID']) in list(unique_RG_IDs):
                if unique_RG_IDs[str(RG_tag['ID'])]['DS'] == str(RG_tag['DS']) and unique_RG_IDs[str(RG_tag['ID'])][
                    'SM'] == str(RG_tag['SM']):
                    RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['outputID'] = str(RG_tag['ID'])
                else:
                    not_unique = True
                    i = 0
                    while not_unique:
                        newID: str = str(RG_tag['ID']) + '_' + str(i)
                        if str(newID) not in list(unique_RG_IDs):
                            RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['outputID'] = str(newID)
                            unique_RG_IDs[str(newID)] = {}
                            unique_RG_IDs[str(newID)]['DS'] = str(RG_tag['DS'])
                            unique_RG_IDs[str(newID)]['SM'] = str(RG_tag['SM'])
                            not_unique = False
                        i += 1
            else:
                RG_ID_dict[str(alignment_file_path)][str(RG_tag['ID'])]['outputID'] = str(RG_tag['ID'])
                unique_RG_IDs[str(RG_tag['ID'])] = {}
                unique_RG_IDs[str(RG_tag['ID'])]['DS'] = str(RG_tag['DS'])
                unique_RG_IDs[str(RG_tag['ID'])]['SM'] = str(RG_tag['SM'])
    return RG_ID_dict, unique_RG_IDs


def _extract_RG_info_from_long_read_input(long_read_input: List[str]) -> object:#tuple[Union[Union[str, bytes], Any], object, Optional[Any], Optional[Any], Optional[Any]]:
    """

    Returns:
        object:
    """
    long_read_input_path = os.path.abspath(long_read_input[0])
    if len(long_read_input) == 1:
        input_ID = None
        input_sample = None
        input_sample_description = None
        input_reference_sequence_input = None
    elif len(long_read_input) == 2:
        input_ID = long_read_input[1]
        input_sample = None
        input_sample_description = None
        input_reference_sequence_input = None
    elif len(long_read_input) == 3:
        input_ID = long_read_input[1]
        input_sample = long_read_input[2]
        input_sample_description = None
        input_reference_sequence_input = None
    elif len(long_read_input) == 4:
        input_ID = long_read_input[1]
        input_sample = long_read_input[2]
        input_sample_description = long_read_input[3]
        input_reference_sequence_input = None
    elif len(long_read_input) >= 5:
        input_ID: object = long_read_input[1]
        input_sample = long_read_input[2]
        input_sample_description = long_read_input[3]
        input_reference_sequence_input = long_read_input[4]
    return long_read_input_path, input_ID, input_sample, input_sample_description, input_reference_sequence_input


def _sample_to_vcf_file_dict(vcf_file_paths):
    """
    Create a dictionary relating sample names to vcf file paths.
    """
    sample_to_vcf_file_dict = {}
    for vcf_file_path in vcf_file_paths:
        vcf_file = pysam.VariantFile(list(vcf_file_path.keys())[0])
        for sample in vcf_file.header.samples:
            sample_to_vcf_file_dict[str(sample)] = vcf_file_path
        vcf_file.close() # Clean up after ourselves!
    return sample_to_vcf_file_dict


def _compile_read_groups(
        alignment_file_path, sample, ID, sample_description, RG_ID_dict, unique_RG_IDs, ignore_samples
        ):
    if ignore_samples:
        RG_ID_dict[str(alignment_file_path)] = 'ignore_samples'
    else:
        # Check integrity of alignment file
        cmd = "samtools quickcheck " + alignment_file_path
        exit_status = os.system(cmd)
        if exit_status != 0:
            # Bam file could not be read.
            raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_compile_read_groups_314)\n" % (alignment_file_path))
            #pass
        
        with pysam.AlignmentFile(alignment_file_path, 'rb') as bam_file:
            RG_tags = bam_file.header.get('RG')
        if sample is not None:
            RG_ID_dict[str(alignment_file_path)] = {}
            if ID is not None:
                RG_ID_dict[str(alignment_file_path)][str(ID)] = {}
                RG_ID_dict[str(alignment_file_path)][str(ID)]['SM'] = str(sample)
                RG_ID_dict[str(alignment_file_path)][str(ID)]['outputID'] = str(ID)
                RG_ID_dict[str(alignment_file_path)][str(ID)]['RG_tags'] = False
                if sample_description is not None:
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(sample_description)
                else:
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(
                            'LRphase_input_file_' + str(alignment_file_path)
                            )
            else:
                ID = '0' + str(random.randint(1, 10000))
                RG_ID_dict[str(alignment_file_path)][str(ID)] = {}
                RG_ID_dict[str(alignment_file_path)][str(ID)]['SM'] = str(sample)
                RG_ID_dict[str(alignment_file_path)][str(ID)]['outputID'] = str(ID)
                RG_ID_dict[str(alignment_file_path)][str(ID)]['RG_tags'] = False
                if sample_description is not None:
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(sample_description)
                else:
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(
                            'LRphase_input_file_' + str(alignment_file_path)
                            )
            if str(ID) not in unique_RG_IDs:
                unique_RG_IDs[str(ID)] = {}
                unique_RG_IDs[str(ID)]['DS'] = str(sample_description)
                unique_RG_IDs[str(ID)]['SM'] = str(sample)
        elif RG_tags is not None:
            RG_ID_dict, unique_RG_IDs = _unique_RG_IDs_from_RG_tags(
                    RG_ID_dict, unique_RG_IDs, alignment_file_path
                    )
        else:
            sys.stderr.write("%s Has No RG tags and was not input with sample information. Reads in this file will not be processed. Either re-input this read file with sample information or resubmit with ignore_samples option.\n" % alignment_file_path)
    return RG_ID_dict, unique_RG_IDs


class InputData:
    """
    sample_to_RG_header: defaultdict[Any, list]
    sample_to_PG_header: defaultdict[Any, list]
    unique_RG_IDs: Dict[Any, Any]
    """
    output_directory: object
    try:
        from LRphase import urls
        urls_found = True
    except:
        sys.stderr.write('Could not find import urls from data. Data will not be able to downloaded from example web sources.\n')
        urls_found = False
    
    def __init__(
            self, output_directory_path: str = None, vcf_file_input: str = None, long_read_input: object = None,
            reference_sequence_input: str = None, sample: str = None, ID: str = None,
            sample_description: str = None,
            ignore_phase_sets: bool = None, ignore_samples: bool = None, download_from_urls: bool = False,
            reference_sequence_input_assembly: str = None, auto_simulate_samples: bool = False,
            max_threads: int = 3, quiet = False, silent = False
            ) -> object:
        """

        Returns:
            object: 
        """
        sys.stderr.write('\n############## Preparing output directory tree ##############\n')
        if output_directory_path is not None:
            self.output_directory = _prepare_output_directory(output_directory_path)
        else:
            self.output_directory = _prepare_output_directory(
                    'LRphase_output_' + str(time.localtime()[0]) + '-' + str(time.localtime()[1]) + '-' + str(
                            time.localtime()[2]
                            ) + '_' + str(time.localtime()[3]) + 'hr_' + str(time.localtime()[4]) + 'min_' + str(
                            time.localtime()[5]
                            ) + 'sec'
                    )
        
        _prepare_output_directory(self.output_directory + '/reference_sequences')
        _prepare_output_directory(self.output_directory + '/haplotype_information')
        _prepare_output_directory(self.output_directory + '/output_reads')
        _prepare_output_directory(self.output_directory + '/input_reads')
        
        #self.summary_file_path = _prepare_summary_file(self.output_directory)
        self.reference_sequence_paths = None
        self.RG_ID_dict = {}
        self.unique_RG_IDs = {}
        self.vcf_files = []
        self.sample_to_vcf_file_dict = {}
        self.sample_to_reference_sequences_dict = defaultdict(list)
        self.sample_to_reference_sequence_path = defaultdict(list)
        self.sample_to_sam_header = defaultdict(list)
        self.sample_to_PG_header = defaultdict(list)
        self.sample_to_RG_header = defaultdict(list)
        self.alignment_file_to_reference_sequences_dict = defaultdict(list)
        self.alignment_files = []
        self.phasable_samples = {}
        self.reference_sequence_files = {}
        
        self.auto_simulate_samples = auto_simulate_samples
        
        # self.sample_hap_to_true_alignment_dict = {}
        
        self.reference_sequence_input = reference_sequence_input
        self.vcf_file_input = vcf_file_input
        self.long_read_input = long_read_input
        self.sample = sample
        self.ID = ID
        self.sample_description = sample_description
        self.ignore_phase_sets = ignore_phase_sets
        self.ignore_samples = ignore_samples
        self.download_from_urls = download_from_urls
        self.reference_sequence_input_assembly = reference_sequence_input_assembly

        self.max_threads = max_threads
        self.quiet = quiet
        self.silent = silent

        # To store log-likelihood ratios and associated error rates from simulations.
        self.lr_err_rate_dict = {}
        
        if not self.reference_sequence_input_assembly is None:
            sys.stderr.write('\n############## Preparing reference sequence directory ##############\n')
            _prepare_output_directory(
                    self.output_directory + '/reference_sequences/' + str(reference_sequence_input_assembly)
                    )

        if not self.long_read_input is None:
            # self.long_read_input = long_read_input
            long_read_input_path, input_ID, input_sample, input_sample_description, input_reference_sequence_input = _extract_RG_info_from_long_read_input(long_read_input)
            if input_reference_sequence_input is None:
                input_reference_sequence_input = self.reference_sequence_input
            self.add_reads(self.long_read_input, input_sample, input_ID,
                           input_sample_description, input_reference_sequence_input)
        
    
    def _download_file(self, url, output_directory = None):
        if output_directory is None:
            local_filename = self.output_directory + '/' + url.split('/')[-1]
        else:
            local_filename = output_directory + '/' + url.split('/')[-1]
        with requests.get(url, stream = True) as r:
            with open(local_filename, 'wb') as f:
                copyfileobj(r.raw, f)
        return local_filename
    

    # Removed commented functions. -- AGD
    

    def add_reference_sequence(
            self, reference_sequence_input, output_directory = None, reference_sequence_input_assembly = None
            ) -> object:
        """
        Create an pysam.FastaFile object and store it within self.reference_sequence_files[].
        First downloads the reference sequence if supplied as a URL.
        If fasta file is bgzipped, unzips it.
        Stores absolute path to the reference sequence file in self.reference_sequence_files[].
        Returns absolute path to current reference sequence file.
        """
        #sys.stderr.write("%s; %s; %s\n" % (reference_sequence_input, output_directory, reference_sequence_input_assembly))

        # First, build the directory path where the reference sequence file will be written.
        if output_directory is None:
            output_directory_path = self.output_directory
        else:
            output_directory_path = output_directory
        
        if reference_sequence_input_assembly is None:
            output_directory_path = self.output_directory + '/reference_sequences'
        else:
            output_directory_path = self.output_directory + '/reference_sequences/' + str(reference_sequence_input_assembly)
            sys.stderr.write('\n############## Preparing reference sequence directory ##############\n')
            _prepare_output_directory(output_directory_path)

        # Download the reference sequence if a URL was provided. Reassign reference_sequence_input
        # to hold the location + name of the downloaded file.
        if reference_sequence_input.startswith('http'):
            sys.stderr.write('The reference sequence input is a url. Downloading now.\n')
            reference_sequence_input = self._download_file(reference_sequence_input, output_directory_path)

        # So this looks like it just reassembles the path to the given input file
        # and applies the .fa extension. Could we just test for different file types
        # instead?
        #reference_sequence_input_path = os.path.dirname(os.path.abspath(reference_sequence_input)) + '/' + os.path.basename(reference_sequence_input).split('.')[0] + '.fa'
        reference_sequence_input_path = '.'.join(reference_sequence_input.split('.')[0:-1]) + '.fa'  # Does same thing with less ops

        # Instantiate the pysam.FastaFile object given the supplied/downloaded file. Assign
        # this object to self.reference_sequence_files, indexed by the file path.
        try:
            ref_seq = pysam.FastaFile(reference_sequence_input_path)
        except (OSError, FileNotFoundError, IsADirectoryError, PermissionError) as e:
            # There was a problem accessing the reference sequence file.
            sys.stderr.write("ERROR: Reference sequence file {} is absent or unreadable.\nError message was: {}\n".format(reference_sequence_input_path, e))
            exit(2)

        # Only attempt to decode if we have a valid, readable file.
        try:
            reference_sequence_file_path: object = ref_seq.filename.decode()
            self.reference_sequence_files[reference_sequence_file_path] = ref_seq
            if reference_sequence_input_assembly:
                self.reference_sequence_files[reference_sequence_input_assembly] = ref_seq

        # If that fails, try unzipping the file with bgzip before instantiating the pysam.FastaFile
        # and making the subsequent value assignments.
        except OSError:
            # This appears to make the implicit assumption that an OSError means
            # that reference_sequence_input is a bgzip file -- maybe not always true??
            # Also, all this does is decompress the input file. Is there a way to
            # avoid this necessity? -- AGD
            sys.stderr.write("ERROR: {}\n".format(OSError))
            with pysam.BGZFile(reference_sequence_input, 'r') as infile:
                with open(reference_sequence_input_path, 'w') as outfile:
                    outfile.write(infile.read().decode())
            ref_seq = pysam.FastaFile(reference_sequence_input_path)
            #reference_sequence_names = ref_seq.references
            reference_sequence_file_path: object = ref_seq.filename.decode()
            self.reference_sequence_files[reference_sequence_file_path] = ref_seq
            if reference_sequence_input_assembly:
                self.reference_sequence_files[reference_sequence_input_assembly] = ref_seq

        # Removed large block of commented code. --AGD

        # Store the reference sequence path to self.reference_sequence_paths[].
        if self.reference_sequence_paths is None:
            self.reference_sequence_paths = [reference_sequence_file_path]
        elif isinstance(self.reference_sequence_paths, str):
            _reference_sequence_paths = []
            _reference_sequence_paths.append(self.reference_sequence_paths)
            self.reference_sequence_paths = _reference_sequence_paths
            self.reference_sequence_paths.append(reference_sequence_file_path)
        elif isinstance(self.reference_sequence_paths, list):
            self.reference_sequence_paths.append(reference_sequence_file_path)
        sys.stderr.write('Reference sequence prepared successfully: %s\n' % reference_sequence_file_path)

        # Return the current reference sequence path.
        return reference_sequence_file_path  #reference_sequence_names

    
    def add_haplotype_information(
            self, vcf_file_input: str, ignore_phase_sets: bool = False, reference_sequence: object = 'hg38'
            ) -> object:
        """
        
        """
        # Check vcf_input_file for proper format and indexing
        vcf_file_path = self._prepare_vcf_file(self.output_directory, vcf_file_input)

        # Check for phasing information
        if pysam.VariantFile(vcf_file_path).header.formats.keys().count('PS') == 0:
            ignore_phase_sets = True
            sys.stderr.write("This VCF file does not have the PS subfield. Phase sets will be ignored and all phased variants on the same chromosome (vcf contig) will be considered to be one contiguous haploblock " % vcf_file_path)

        # Append the VCF file path to self.vcf_files[], etc.
        self.vcf_files.append({vcf_file_path:ignore_phase_sets})
        self.sample_to_vcf_file_dict = _sample_to_vcf_file_dict(self.vcf_files)
        self.sample_to_alignment_files = _sample_to_alignment_files(self.sample_to_vcf_file_dict, self.RG_ID_dict)
        try:
            self._sample_to_reference_sequences_dict()
        except Exception as e:
            sys.stderr.write("%s\n" % (e))
            exit(2)
        self._sample_to_PG_dict()

        # Return path to the vcf file.
        return vcf_file_path

    
    def __iter__(self):
        self.sample_counter = 0
        self.phasable_samples = {}
        return self
    
    def new_PhasableSample(self, sample: str, reference_sequence_paths: list = None, auto_simulate_samples: bool = False) -> object:
        """

        Args:
            sample:
            reference_sequence_paths:

        Returns:

        """
        vcf_file_path, ignore_phase_sets = _pair_sample_with_vcf(
                sample, self.sample_to_vcf_file_dict, self.ignore_phase_sets
                )
        reference_sequence_names: list = self.sample_to_reference_sequences_dict[sample]
        #sys.stderr.write("%s\n" % sample)
        #sys.stderr.write("%s\n" % self.sample_to_alignment_files[sample])
        phasable_sample = PhasableSample.PhasableSample(
            sample, vcf_file_path, ignore_phase_sets, self.sample_to_alignment_files[sample], self.RG_ID_dict,
            reference_sequence_names, reference_sequence_paths = reference_sequence_paths,
            simulate_haplotypes = auto_simulate_samples, output_directory_path=self.output_directory,
            quiet = self.quiet, silent = self.silent
        )
        return phasable_sample
    
    def __next__(self):
        if self.sample_counter < len(list(self.sample_to_vcf_file_dict.keys())):
            sample = str(list(self.sample_to_vcf_file_dict.keys())[self.sample_counter])
            if sample in self.sample_to_reference_sequence_path:
                if len(self.sample_to_reference_sequence_path[sample]) > 0:
                    reference_sequence_paths = self.sample_to_reference_sequence_path[sample]
            else:
                reference_sequence_paths = self.reference_sequence_paths
            phasable_sample = self.new_PhasableSample(sample, reference_sequence_paths, self.auto_simulate_samples)
            #vcf_file_path, ignore_phase_sets = self._pair_sample_with_vcf(sample, self.sample_to_vcf_file_dict,
            #                                                             self.ignore_phase_sets)
            #phasable_sample = PhasableSample(sample, vcf_file_path, ignore_phase_sets,
            #                               self.sample_to_alignment_files[sample], self.RG_ID_dict)
            self.sample_counter += 1
            self.phasable_samples[phasable_sample.sample] = phasable_sample
            return phasable_sample
        else:
            raise StopIteration()


    """
    As implemented, this holds all the reads in memory simultaneously. As of 12/15/2022, this leads to 
    excessive memory usage (>100G) and is causing problems for users. It would be better to structure this
    so that it reads-in/holds in memory only what it needs at any given time.
    """
    def add_reads(
            self, long_reads_alignment_path, sample = None, ID = None, haplotype = None, sample_description = None,
            reference_sequence_input = None, database = None, master_database = False, simulated = False,
            reference_sequence_input_assembly = None, preserve_alignment = False, clobber = False,
            quiet = False, silent = False, threads = 1, no_align = False
            ):
        """
        clobber arg will cause any existing reads/alignments to be purged prior to loading current read set when True.
        """
        if clobber:
            self._reinit_reads(preserve_alignments=preserve_alignment)
        
        if long_reads_alignment_path.startswith('http'):
            sys.stderr.write('The reads input is a url. Downloading now.\n')
            long_reads_alignment_path = self._download_file(long_reads_alignment_path, self.output_directory)

        try:
            sorted_bam_file_paths, combined_long_read_fastq_path = self._parse_long_reads_input(long_reads_alignment_path,
                                                                                                self.output_directory,
                                                                                                threads=threads)
        except Exception as e:
            sys.stderr.write("Error parsing long-read file: %s\n" % e)
            exit(2)
        if combined_long_read_fastq_path:
            if reference_sequence_input is None and self.reference_sequence_input is None:
                sys.stderr.write("Must supply reference sequence!\n")
                raise Exception("Reference sequence not given.")
            elif self.reference_sequence_input is None:
                reference_sequence_input = self.add_reference_sequence(
                    reference_sequence_input = reference_sequence_input,
                    output_directory = self.output_directory,
                    reference_sequence_input_assembly = reference_sequence_input_assembly
                )
            else:
                reference_sequence_input = self.reference_sequence_input
                
            sorted_bam_file_paths.append(
                _align_long_reads_fastq(combined_long_read_fastq_path,
                                        reference_sequence_input,
                                        self.output_directory,
                                        self.max_threads,
                                        quiet = quiet,
                                        silent = silent,
                                        no_align = no_align
                )
            )
        for _sorted_bam_file_path in sorted_bam_file_paths:
            if _sorted_bam_file_path:
                sorted_bam_file_path, combined_long_read_fastq_path = self._parse_long_reads_input(
                    _sorted_bam_file_path,
                    self.output_directory,
                    threads = threads,
                    is_sorted_bam = no_align
                )
                #sys.stderr.write("%s; %s\n" % (sorted_bam_file_path, combined_long_read_fastq_path))
                self.RG_ID_dict, self.unique_RG_IDs = _compile_read_groups(
                    sorted_bam_file_path[0], sample, ID, sample_description,
                    self.RG_ID_dict, self.unique_RG_IDs,
                    self.ignore_samples
                )
                #sys.stderr.write("%s; %s\n" % (self.RG_ID_dict, self.unique_RG_IDs))
                self.alignment_files.append(sorted_bam_file_path[0])
                #sys.stderr.write("%s\n" % self.alignment_files)
                self.sample_to_alignment_files = _sample_to_alignment_files(
                        self.sample_to_vcf_file_dict, self.RG_ID_dict
                        )
                #sys.stderr.write("%s\n" % self.sample_to_alignment_files)
                self._sample_to_reference_sequences_dict(sorted_bam_file_path[0])
                self._sample_to_PG_dict(sam_file_path = sorted_bam_file_path[0])
                if reference_sequence_input is not None:
                    if not reference_sequence_input in self.sample_to_reference_sequence_path[sample]:
                        self.sample_to_reference_sequence_path[sample].append(reference_sequence_input)


    def _reinit_reads(self, preserve_alignments=True):
        """
        Reinitialize all objects related to storing reads in an InputData objects.
        """
        self.RG_ID_dict = {}
        self.unique_RG_IDs = {}
        if not preserve_alignments:
            self._purge_alignment_files()
        self.alignment_files = []
        self.sample_to_reference_sequences_dict = defaultdict(list)
        self.sample_to_reference_sequence_path = defaultdict(list)
        self.sample_to_PG_header = defaultdict(list)
        self.sample_to_RG_header = defaultdict(list)
        return
    

    def _purge_alignment_files(self):
        """
        Purge alignment files from the filesystem.
        """
        for alignment_file in self.alignment_files:
            if os.path.isfile(alignment_file):
                os.remove(alignment_file)
            if os.path.exists("%s.bai" % (alignment_file)):
                os.remove("%s.bai" % (alignment_file))            
        self.alignment_files = []
        return

    
    def _prepare_sam_header(self, _sample = None):
        if not str(_sample) in self.sample_to_sam_header:
            self.sample_to_sam_header[_sample] = defaultdict(list)
        self.sample_to_sam_header[_sample]['HD'] = {'VN':'1.6', 'SO':'coordinate'}
        self.sample_to_sam_header[_sample]['SQ'] = self.sample_to_reference_sequences_dict[_sample]
        self.sample_to_sam_header[_sample]['RG'] = self.sample_to_RG_header[_sample]
        self.sample_to_sam_header[_sample]['PG'] = self.sample_to_PG_header[_sample]
        return  # self.sample_to_sam_header[_sample]
    

    def _sample_to_reference_sequences_dict(self, sam_file_path = None, sam_file = None):
        """
        Prepare a dictionary relating reference sequences to samples.
        """
        #if sam_file_path == None and sam_file == None:
        #    raise Exception("No valid alignment files given. IP_sample_to_ref_seq_dict_753")

        if sam_file_path:
            # Check integrity of alignment file
            cmd = "samtools quickcheck " + sam_file_path
            exit_status = os.system(cmd)
            if exit_status != 0:
                # Bam file could not be read.
                raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_sample_to_refseq_dict_764)\n" % (sam_file_path))
                #pass
        
            aln_file = pysam.AlignmentFile(sam_file_path)
            for reference_sequence in aln_file.header['SQ']:
                self.alignment_file_to_reference_sequences_dict[sam_file_path].append(reference_sequence)
            aln_file.close()
        for sample in self.sample_to_vcf_file_dict:
            for sam_file_path in self.alignment_file_to_reference_sequences_dict:
                if sam_file_path in self.sample_to_alignment_files[str(sample)]:
                    # Check file integrity
                    cmd = "samtools quickcheck " + sam_file_path
                    exit_status = os.system(cmd)
                    if exit_status != 0:
                        # Bam file could not be read.
                        sys.stderr.write("WARNING: %s could not be read. Reads will not be analysed. This is likely due to a file that is incomplete or corrupt.\n" % (sam_file_path))
                        continue
                    
                    aln_file = pysam.AlignmentFile(sam_file_path)
                    for reference_sequence in aln_file.header['SQ']:
                        if str(reference_sequence['SN']) not in [str(reference_seq['SN']) for reference_seq in
                                                                 self.sample_to_reference_sequences_dict[sample]]:
                            self.sample_to_reference_sequences_dict[sample].append(reference_sequence)
                    aln_file.close()
                    
    
    def _sample_to_PG_dict(self, sam_file_path: object = None, sam_file: object = None) -> object:
        """
        Creates a dictionary relating samples to phasing groups.
        Args:
            sam_file_path:
            sam_file:
        """
        #if sam_file_path == None and sam_file == None:
        #    raise Exception("No valid alignment files given. IP_sample_to_PG_dict_798")
        
        if sam_file_path:
            # Check file integrity
            cmd = "samtools quickcheck " + sam_file_path
            exit_status = os.system(cmd)
            if exit_status != 0:
                # Bam file could not be read.
                raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_sample_to_PG_dict_805)\n" % (sam_file_path))
                #pass
            
            aln_file = pysam.AlignmentFile(sam_file_path)
            for reference_sequence in aln_file.header['SQ']:
                self.alignment_file_to_reference_sequences_dict[sam_file_path].append(reference_sequence)
            aln_file.close()
            
        for sample in self.sample_to_vcf_file_dict:            
            for sam_file_path in self.alignment_file_to_reference_sequences_dict:
                if sam_file_path in self.sample_to_alignment_files[str(sample)]:
                    # Check file integrity
                    cmd = "samtools quickcheck " + sam_file_path
                    exit_status = os.system(cmd)
                    if exit_status != 0:
                        # Bam file could not be read.
                        sys.stderr.write("WARNING: %s could not be read. Reads will not be analysed. This is likely due to a file that is incomplete or corrupt.\n" % (sam_file_path))
                        continue

                    aln_file = pysam.AlignmentFile(sam_file_path)
                    for PG_tag in aln_file.header['PG']:
                        PG_tag['ID'] = str(list(self.RG_ID_dict[sam_file_path].keys())[0])
                        # if PG_tag['ID'] in [PG_tag['ID'] for PG_tag in self.sample_to_PG_header[sample]]
                        RG_tag = {
                                'ID':str(list(self.RG_ID_dict[sam_file_path].keys())[0]), 'SM':str(
                                        self.RG_ID_dict[sam_file_path][
                                            str(list(self.RG_ID_dict[sam_file_path].keys())[0])]['SM']
                                        ),
                                'DS':str(
                                        self.RG_ID_dict[sam_file_path][
                                            str(list(self.RG_ID_dict[sam_file_path].keys())[0])]['DS']
                                        )
                                }
                        # {str(list(self.RG_ID_dict[sam_file_path].keys())[0]): {'DS': 'small7.fastq', 'SM': 'HG001',
                        # 'RG_tags': True, 'outputID': '1'}} self.RG_ID_dict[sam_file_path] self.RG_ID_dict[
                        # sam_file_path][str(list(self.RG_ID_dict[sam_file_path].keys())[0])]['SM'] if not PG_tag[
                        # 'ID'] in [PG_tag['ID'] for PG_tag in self.sample_to_PG_header[sample]]:
                        # self.sample_to_PG_header[sample].append(PG_tag)
                        new_PG_tag: bool = True
                        for PG_header in self.sample_to_PG_header[sample]:
                            if PG_tag['PN'] == PG_header['PN'] and PG_tag['ID'] == PG_header['ID']:
                                new_PG_tag = False
                        if new_PG_tag:
                            self.sample_to_PG_header[sample].append(PG_tag)
                        # self.sample_to_PG_header[sample].append(PG_tag)
                        if not RG_tag['ID'] in [RG_tag['ID'] for RG_tag in self.sample_to_RG_header[sample]]:
                            self.sample_to_RG_header[sample].append(RG_tag)
                    aln_file.close()

    
    def _parse_long_reads_input(self, long_read_input: object, output_directory: str, threads: int = 1, is_sorted_bam = False) -> object:
        """

        Args:
            long_read_input (object):
            output_directory:

        Returns:
            object:

        """
        if long_read_input.startswith('http'):
            sys.stderr.write('The long read input is a url. Downloading now.\n')
            long_reads_alignment_path = self._download_file(long_read_input, self.output_directory + '/input_reads')
        else:
            long_reads_alignment_path = long_read_input
        
        combined_long_read_fastq_path = [] # Why is this an array? It gets assigned a string value by os.path.abspath.
        sorted_bam_file_paths = []

        if is_sorted_bam:
            sorted_bam_file_path = os.path.abspath(
                '%s/%s%s' % (
                    output_directory,
                    os.path.splitext(os.path.basename(long_reads_alignment_path))[0],
                    '_sorted.bam'
                )
            )
            sorted_bam_file_paths.append(sorted_bam_file_path)
            return sorted_bam_file_paths, combined_long_read_fastq_path
            
        elif not os.path.exists(long_reads_alignment_path):
            sys.stderr.write("ERROR: Could not find %s.\nUse -i to specify the path of a file containing reads for phasing OR use -i to specify the path of a directory containing the long read files and all files will be processed individually.\n(EX: -i path/minion_run3_GM12878/minion_run3_GM12878_0.fastq OR -i path/minion_run3_GM12878/minion_run3_GM12878_0.sam OR -i path/minion_run3_GM12878/minion_run3_GM12878_0_sorted.bam OR -i path/minion_run3_GM12878/\n\n)." % long_reads_alignment_path)
            raise OSError("File not found")
        
        elif os.path.isdir(long_reads_alignment_path):
            sys.stderr.write('Directory was given as input for read files. Processing all files in %s.\n' % long_reads_alignment_path)
            for aln_file in os.listdir(long_reads_alignment_path):
                qualified_aln_file = '%s/%s' % (os.path.abspath(long_reads_alignment_path), aln_file)
                if os.path.isfile(qualified_aln_file):
                    if str(qualified_aln_file).lower().endswith('.fastq') or str(qualified_aln_file).lower().endswith('.fastq.gz'):
                        sys.stderr.write('%s is a valid fastq file.\n' % qualified_aln_file)
                        if not combined_long_read_fastq_path:
                            combined_long_read_fastq_path = os.path.abspath(
                                    '%s/%s%s' % (
                                            output_directory,
                                            os.path.splitext(os.path.basename(long_reads_alignment_path))[0],
                                            '_combined_fastq.gz')
                                    )
                        with open(combined_long_read_fastq_path, 'a') as combined_fastqs:
                            with open(qualified_aln_file, 'r') as fastqfile:
                                combined_fastqs.write(fastqfile.read())
                    
                    elif str(qualified_aln_file).lower().endswith('.sam') or str(file).lower().endswith('.bam'):
                        sorted_bam_file_paths.append(
                            _prepare_alignment(
                                output_directory,
                                file,
                                threads=threads                            )
                        )
        
        elif os.path.isfile(long_reads_alignment_path):
            if str(long_reads_alignment_path).lower().endswith('.fastq') or str(
                    long_reads_alignment_path
                    ).lower().endswith('.fastq.gz'):
                sys.stderr.write('%s is a valid fastq file.\n' % long_reads_alignment_path)
                combined_long_read_fastq_path = os.path.abspath(long_reads_alignment_path)
            
            elif str(long_reads_alignment_path).lower().endswith('.sam') or str(
                    long_reads_alignment_path
                    ).lower().endswith('.bam'):
                sorted_bam_file_paths.append(
                    _prepare_alignment(
                        output_directory,
                        long_reads_alignment_path,
                        threads=threads
                    )
                )
        else:
            sys.stderr.write('Error: Reads should be in .fastq, .fastq.gz, .sam, or .bam format. %s does not have a correct suffix to be a valid format and is not a directory. Use -i to specify the path of a file containing reads for phasing OR use -i to specify the path of a directory containing the long read files and all files will be processed individually. (EX: -i path/minion_run3_GM12878/minion_run3_GM12878_0.fastq OR -i path/minion_run3_GM12878/minion_run3_GM12878_0.sam OR -i path/minion_run3_GM12878/minion_run3_GM12878_0_sorted.bam OR -i path/minion_run3_GM12878/).' % long_reads_alignment_path)
            raise Exception('Unsupported long read file format.')
        
        return sorted_bam_file_paths, combined_long_read_fastq_path
    
    def _prepare_vcf_file(self, output_directory: object, vcf_file_input: str) -> object:
        """
        Checks the variant file for proper GT headers, bgzip format, and tabix indexing.
        Attempts to sort/bgzip/index if needed.
        Returns path to bgzipped/indexed vcf.
        """
        start_process_time = time.time()
        sys.stderr.write('\n############## Preparing vcf file ##############\n')
        
        if vcf_file_input.startswith('http'):
            sys.stderr.write('The vcf file input is a url. Downloading now.\n')
            vcf_file_input = self._download_file(vcf_file_input, self.output_directory + '/haplotype_information')
        
        if not os.path.isfile(vcf_file_input):
            sys.stderr.write('Could not find %s. Use -v to specify the path of the vcf file to be used as haplotype information for phasing. (EX: -v path/to/GM12878.vcf.gz or --vcf GM12878.vcf).' % vcf_file_input)
            return

        # Check input file for proper format (GT header present, bgzipped, tabix indexed).
        # Correct the file format and index if needed.
        vcf_file_pysam: VariantFile = pysam.VariantFile(vcf_file_input)
        if vcf_file_pysam.format == 'VCF': # What if it's not???
             # Check to be sure the GT header is present.
            if vcf_file_pysam.header.formats.keys().count('GT') == 0:
                sys.stderr.write('The VCF file provided does not have the GT subfield. VCF files must have the GT subfield for all samples in order to extract genotype information. Phased variants should have | in the GT subfield instead of / \n')
                return
            
            # Check for bgzip compression.
            if vcf_file_pysam.compression == 'BGZF':
                vcf_file_path = os.path.abspath(vcf_file_input)
                sys.stderr.write('%s is a valid vcf file.\n' % vcf_file_path)
            # Compress with bgzip if not already done.
            else:
                sys.stderr.write('%s is a valid .vcf file but it is not in bgzip (.vcf.gz) format. VCF files must be compressed with bgzip and indexed with tabix. LRphase will now attempt to run bgzip on %s.\n' % (vcf_file_input, vcf_file_input))
                # Assuming there should be a call to _sort_vcf_file here??
                vcf_file_path = _sort_vcf_file(vcf_file_path)

            # Check for a tabix index and create one if not found.
            vcf_file_index_path = vcf_file_path + '.tbi'    
            if os.path.isfile(vcf_file_index_path):
                sys.stderr.write('Found %s as an index for vcf file' % vcf_file_index_path)
            else:
                sys.stderr.write('%s is a valid .vcf file in bgzip (.vcf.gz) format but an index was not found. Indexing with tabix now.\n' % vcf_file_input)
                pysam.tabix_index(vcf_file_path, preset = 'vcf', force = True)
            # The VariatnFile object is not stored for later use, so close it to free up resources.
            vcf_file_pysam.close()
        else:
            sys.stderr.write("ERROR: %s is not a valid VCF file.\n" % vcf_file_input)
            return

        # Calculate and report time elapsed.
        total_process_time = time.time() - start_process_time       
        sys.stderr.write('Prepared vcf file in %.2f seconds\n' % total_process_time)

        # Return the path to the bgzipped and indexed vcf.
        return vcf_file_path

    
    def _compile_read_groups(
            self, alignment_file_path, sample, ID, sample_description, RG_ID_dict, unique_RG_IDs, ignore_samples
    ):
        if ignore_samples:
            RG_ID_dict[str(alignment_file_path)] = 'ignore_samples'
        else:
            # Check file integrity before opening
            cmd = "samtools quickcheck " + alignment_file_path
            exit_status = os.system(cmd)
            if exit_status != 0:
                # Bam file could not be read.
                raise Exception("\nERROR: %s could not be read. This is likely due to a file that is incomplete or corrupt. Exiting. (IDP_compile_read_groups_1007)\n" % (alignment_file_path))
                #pass
            
            with pysam.AlignmentFile(alignment_file_path, 'rb') as bam_file:
                RG_tags = bam_file.header.get('RG')
            if sample is not None:
                RG_ID_dict[str(alignment_file_path)] = {}
                if ID is not None:
                    RG_ID_dict[str(alignment_file_path)][str(ID)] = {}
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['SM'] = str(sample)
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['outputID'] = str(ID)
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['RG_tags'] = False
                    if sample_description is not None:
                        RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(sample_description)
                    else:
                        RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(
                            'LRphase_input_file_' + str(alignment_file_path)
                        )
                else:
                    ID = '0' + str(random.randint(1, 10000))
                    RG_ID_dict[str(alignment_file_path)][str(ID)] = {}
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['SM'] = str(sample)
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['outputID'] = str(ID)
                    RG_ID_dict[str(alignment_file_path)][str(ID)]['RG_tags'] = False
                    if sample_description is not None:
                        RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(sample_description)
                    else:
                        RG_ID_dict[str(alignment_file_path)][str(ID)]['DS'] = str(
                            'LRphase_input_file_' + str(alignment_file_path)
                        )
                if str(ID) not in unique_RG_IDs:
                    unique_RG_IDs[str(ID)] = {}
                    unique_RG_IDs[str(ID)]['DS'] = str(sample_description)
                    unique_RG_IDs[str(ID)]['SM'] = str(sample)
            elif RG_tags is not None:
                RG_ID_dict, unique_RG_IDs = _unique_RG_IDs_from_RG_tags(
                    RG_ID_dict, unique_RG_IDs, alignment_file_path
                )
            else:
                sys.stderr.write('%s has No RG tags and was not input with sample information. Reads in this file will not be processed. Either re-input this read file with sample information or resubmit with ignore_samples option.\n' % alignment_file_path)
        return RG_ID_dict, unique_RG_IDs


    def _tabulate_sim_results(self):
        """
        Tabulate log-likelihoods and associated error rates for simulation results.

        """
        for phasable_sample in self:
            if phasable_sample.sample == 'HG001':
                for alignment in phasable_sample:
                    phased_read = PhasedRead(alignment, vcf_file = phasable_sample.vcf_file_path, sample = 'HG001', evaluate_true_alignment = True)
                    if phased_read.is_phased_correctly == False:
                        #print(str(phased_read))	
                        #print(str(phased_read.matches_true_alignment))
                        #print("{}\t{}\t{}\n".format(str(phased_read.phase), str(phased_read.is_phased_correctly), str(phased_read.log_likelihood_ratio)))
                        lr = str(phased_read.log_likelihood_ratio)                        
                        if lr in self.lr_err_rate_dict.keys():
                            self.lr_err_rate_dict[lr] += 1
                        else:
                            self.lr_err_rate_dict[lr] = 1
        return

        

def _sort_vcf_file(vcf_file_input: object, vcf_file_output: object = None) -> object:
    """

    Args:
        vcf_file_input: Location of VCF file (absolute path)
        vcf_file_output: Location/name of VCF.bgz output (absolute path)

    Returns path to sorted/bgzipped output file.
    """
    if vcf_file_output == None:
        vcf_file_output = vcf_file_input + '.bgz'
    
    try:
        subprocess.run(
            ["grep -v ^'#' %s | sort -k1,1 -k2,2n | bgzip > %s" % (vcf_file_input, vcf_file_output)],
            check = True, shell = True
        )
        return vcf_file_output
    except Exception as e:
        sys.stderr.write("Error occurred when bgzip was run. Bgzip must be installed on PATH for LRphase to continue. Please check bgzip installation or provide a bgzip compressed vcf file using the -v option (EX: -v path/to/GM12878.vcf.gz). Error: %s\n" % e)
        return

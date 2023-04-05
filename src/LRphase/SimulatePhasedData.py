# coding=utf-8
import sys, os, re
import subprocess
import time
from glob import glob
import pyliftover
from pysam import faidx, bcftools, FastxFile, AlignedSegment
from typing import Tuple, List


def simulate_reads_pbsim2(
        reference_sequence: str, path_to_pbsim: str = 'pbsim', depth: int = 1,
        simulation_mode: str = 'pbsim2/data/R103.model',
        difference_ratio: str = '23:31:46', length_mean: int = 20000,
        length_max: int = 1000000, length_min: int = 100,
        length_sd: int = 15000, accuracy_min: float = 0.01,
        accuracy_max: float = 1.00, accuracy_mean: float = 0.98,
        prefix: str = None, id_prefix: str = 'S', output_directory: str = None, sample: str = None,
        haplotype: int = None, quiet = False, silent = False, no_sim = False
        ) -> str:
    """
    Args:
        reference_sequence (str):
            reference sequence template to simulate reads from
        path_to_pbsim (str):
            if pbsim2 is not installed to the PATH as pbsim then use this to set its absolute or relative path
        depth (int):
            the depth of coverage that will be simulated
        simulation_mode: 
        difference_ratio: 
        length_mean: 
        length_min: 
        length_max: 
        length_sd: 
        accuracy_min: 
        accuracy_max: 
        accuracy_mean: 
        prefix: 
        id_prefix: 
        output_directory: 
        sample: 
        haplotype: 

    Returns:
        object: 
    
    """
    if output_directory is None:
        output_directory: str = 'simulated'
    
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    
    if not sample is None and not haplotype is None:
        if prefix is None:
            prefix = str(output_directory) + '/' + str(sample) + '__' + str(haplotype) + '__' + 'pbsim2_simulated_' +\
                     os.path.basename(reference_sequence).split('.')[0] + '_' + str(time.localtime()[0]) + '_' + str(
                    time.localtime()[1]
                    ) + '_' + str(time.localtime()[2]) + '_' + str(time.localtime()[3]) + 'hr_' + str(
                    time.localtime()[4]
                    ) + 'min_' + str(time.localtime()[5]) + 'sec'
        else:
            prefix = '{0}/{1}__{2}__{3}'.format(str(output_directory), str(sample), str(haplotype), str(prefix))
    else:
        prefix = '{0}/{1}__{2}__{3}'.format(str(output_directory), str(sample), str(haplotype), str(prefix))
    
    if str(simulation_mode).lower().endswith(".fastq"):
        mode1 = '--sample-fastq'
        mode2 = simulation_mode
    
    elif str(simulation_mode).lower().endswith('.model'):
        mode1 = '--hmm_model'
        mode2 = simulation_mode
    
    if not reference_sequence:
        sys.stderr.write('No reference genome was provided for pbsim2 simulation. Please specify reference_genome_path.\n')
        raise exception('Reference sequence not provided.')
        return
    
    sys.stderr.write('\n################ Beginning Read Simulation ################\n')
    final_combined_simulated_long_reads_fastq_path = '{0}/{1}.fastq'.format(
        str(output_directory),
        str(os.path.splitext(os.path.basename(reference_sequence))[0])
    )
    if no_sim:
        return final_combined_simulated_long_reads_fastq_path
    
    start_process_time = time.time()
    
    reference_genome_path = os.path.abspath(reference_sequence)
    reference_genome_index_path = str(reference_genome_path) + '.fai'
    if not os.path.isfile(reference_genome_index_path):
        faidx(reference_genome_path)

    pbsim2_args = [str(path_to_pbsim), f'--depth', str(depth),
                   str(mode1), str(mode2),
                   f'--difference-ratio', str(difference_ratio),
                   f'--length-mean', str(length_mean),
                   f'--length-min', str(length_min),
                   f'--length-max', str(length_max),
                   f'--length-sd', str(length_sd),
                   f'--accuracy-mean', str(accuracy_mean),
                   f'--accuracy-min', str(accuracy_min),
                   '--accuracy-max', str(accuracy_max),
                   f'--prefix', str(prefix),
                   '--id-prefix', str(id_prefix),
                   str(reference_genome_path)]

    sys.stderr.write("pbsim2 command: %s\n" % (' ' .join(pbsim2_args)))
    
    if silent:
        subprocess.run(pbsim2_args, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    elif quiet:
        subprocess.run(pbsim2_args, check = True, stdout = subprocess.DEVNULL)
    else:
        subprocess.run(pbsim2_args, check = True)

    sys.stderr.write('Post-processing simulated fastq files...\n')
    sample_hap_to_true_alignment_dict = build_true_alignment_dict(reference_genome_index_path, prefix)

    for fastq_file in glob(str(prefix) + '*.fastq'):
        if not (quiet or silent):
            sys.stderr.write("%s\n" % fastq_file)
        with FastxFile(fastq_file, 'r') as infile:
            with open(final_combined_simulated_long_reads_fastq_path, 'a') as outfile:
                for read in infile:
                    read.name = str(sample) + '__' + str(haplotype) + '__' + str(sample_hap_to_true_alignment_dict[read.name])
                    outfile.write(str(read) + '\n')
        os.remove(fastq_file)

    # Clean up pbsim2 output files
    for f in glob(str(prefix) + '*.maf'):
        os.remove(f)
    for f in glob(str(prefix) + '*.ref'):
        os.remove(f)

    sys.stderr.write('Simulation finished in %.2f seconds.\n' % (time.time() - start_process_time))
    return final_combined_simulated_long_reads_fastq_path


def build_true_alignment_dict(reference_genome_index_path, prefix):
    """
    Build the true alignment dict from pbsim2 maf files.
    """
    true_alignment_dict = {}
    files = glob(prefix + '*.maf')
    for pbsim2_maf in files:
        al_dict = extract_true_alignments_from_pbsim2_maf(reference_genome_index_path, pbsim2_maf)
        true_alignment_dict.update(al_dict)
    return true_alignment_dict


def extract_true_alignments_from_pbsim2_maf(reference_genome_index_path, pbsim2_maf):
    """
    Parse true alignment coordinates from a pbsim2 MAF file.
    """
    # Build the chromosome list from the reference genome index
    chr_list = []
    with open(reference_genome_index_path, 'r') as index_file:
        for rec in index_file:
            line = rec.strip('\n').split()
            if len(line) > 0:
                chr_list.append(line[0])

    # Parse MAF records, extracting the true alignment coordinates for each simulated read.
    true_alignment_dict = {}
    state = 0
    reg = []
    with open(pbsim2_maf) as maf_file:
        for rec in maf_file:
            line = rec.strip('\n').split()
            if len(line) == 0:
                continue
            elif state == 0 and line[0] == 'a':
                state = 1
            elif state == 1 and line[0] == 's':
                reg = [line[2], str(int(line[2])+int(line[3]))]
                state = 2
            elif state == 2 and line[0] == 's':
                m = re.match('^S\d+_\d+', line[1])
                if not m:
                    sys.stderr.write("Failed to parse read name from record: %s. Skipping...\n" % " ".join(line[0:5]))
                    state = 0
                    continue
                chr_id = int(line[1].split('S')[1].split('_')[0]) - 1
                if chr_id > len(chr_list):
                    sys.stderr.write("Chromosome not in index: %s. Skipping...\n" % str(chr_id+1))
                    state = 0
                    continue
                name = "!".join([line[1], chr_list[chr_id], reg[0], reg[1], line[4]])
                seq = line[6].replace('-', '')
                if len(seq) != int(line[5]):
                    sys.stderr.write("Inconsistent read length for sequence %s. Skipping...\n" % name)
                    state = 0
                    continue
                true_alignment_dict[line[1]] = name
                #sys.stderr.write("%s, %s\n" % (line[1], name))
                state = 0

    return true_alignment_dict
    

def convert_coordinates(
        haplotype: int, sample: str, contig: str, position: int, reference_sequence: str,
        reference_liftover_converters = None, chain_file = None
        ):
    """

    Args:
        haplotype (int):
        sample (str):
        contig (str):
        position (int):
        reference_sequence (str):
        reference_liftover_converters: 
        chain_file: 

    Returns:
        object: 

    """
    if reference_liftover_converters is None:
        reference_liftover_converters = {}
    
    if reference_sequence in reference_liftover_converters:
        if sample in reference_liftover_converters[reference_sequence]:
            if str(haplotype) in reference_liftover_converters[reference_sequence][sample]:
                return reference_liftover_converters[reference_sequence][sample][str(haplotype)].convert_coordinates(
                        contig, position
                        )
        
        elif chain_file is not None:
            if os.path.isfile(chain_file):
                reference_liftover_converters[reference_sequence][sample] = {}
                reference_liftover_converters[reference_sequence][sample][str(haplotype)] = pyliftover.LiftOver(
                        chain_file
                        )
                return reference_liftover_converters[reference_sequence][sample][str(haplotype)].convert_coordinates(
                        contig, position
                        )
    
    elif chain_file is not None:
        if os.path.isfile(chain_file):
            reference_liftover_converters[reference_sequence] = {}
            reference_liftover_converters[reference_sequence][sample] = {}
            reference_liftover_converters[reference_sequence][sample][str(haplotype)] = pyliftover.LiftOver(chain_file)
            return reference_liftover_converters[reference_sequence][sample][str(haplotype)].convert_coordinates(
                    contig, position
                    )


def generate_haplotype_specific_fasta(
        haplotype: int, sample: str, input_reference_sequence_path: str, haplotype_vcf: str,
        output_reference_sequence_path: str = None,
        chain_file_path: str = None,
        clobber = False
        ) -> object:
    """

    Args:
        haplotype:
        sample:
        input_reference_sequence_path:
        haplotype_vcf:
        output_reference_sequence_path:
        chain_file_path:

    Returns:

    """
    sys.stderr.write('\n################ Generating Reference Sequence for Haplotype %d ################\n' % haplotype)
    input_reference_sequence_path = os.path.abspath(input_reference_sequence_path)
    input_reference_sequence_index_path: str = str(input_reference_sequence_path) + '.fai'

    if output_reference_sequence_path is None:
        output_reference_sequence_path: str = str(os.path.splitext(input_reference_sequence_path)[0]) + '_hap' + str(
                haplotype
                ) + '_' + str(sample) + '.fa'
    
    if chain_file_path is None:
        chain_file_path = str(os.path.splitext(input_reference_sequence_path)[0]) + '_hap' + str(haplotype) + '_' + str(
                sample
                ) + '.chain'

    if os.path.exists(output_reference_sequence_path) and not clobber:
        sys.stderr.write("%s already exists. Skipping.\n" % (output_reference_sequence_path))
    else:
        # This causes python to segfault if bcftools encounters an error!
        fasta_output = bcftools.consensus(
            '-H', str(haplotype), '-s', str(sample), '-c', str(chain_file_path), '-f',
            str(input_reference_sequence_path),
            str(haplotype_vcf)
        )
        with open(output_reference_sequence_path, 'w') as output_fasta:
            output_fasta.write(fasta_output)

    sys.stderr.write("Done.\n")
    return output_reference_sequence_path, chain_file_path


def index_fasta_file(fasta_file_path: str) -> str:
    """

    Args:
        fasta_file_path: 

    Returns:

    """
    fasta_file_path = os.path.abspath(fasta_file_path)
    
    if fasta_file_path.endswith('.fa') or fasta_file_path.endswith('fa.gz'):
        faidx(fasta_file_path)
        return os.path.abspath(fasta_file_path)
    
    elif fasta_file_path.endswith('fa.gz'):
        os.rename(fasta_file_path, '{0}.fa.gz'.format(str(os.splitext(fasta_file_path)[0])))
        faidx(str(os.splitext(fasta_file_path)[0]) + '.fa.gz')
        return str(os.splitext(fasta_file_path)[0]) + '.fa.gz'
    
    else:
        os.rename(fasta_file_path, str(os.splitext(fasta_file_path)[0]) + '.fa')
        faidx(str(os.splitext(fasta_file_path)[0]) + '.fa')
        return str(os.splitext(fasta_file_path)[0]) + '.fa'


def create_fasta_file(
        input_fasta_file_path: str, output_fasta_file_path: str = None, regions: List[str] = None, only_autosomal:
        bool = False, chrom_sizes: str = None
        ) -> str:
    """

    Args:
        input_fasta_file_path: 
        only_autosomal (object): 
        regions:
            List of regions and/or contigs ie: ['chr1:10000000:60000000',...]
        output_fasta_file_path (object): 
    """
    regions_file: str = str(os.path.splitext(output_fasta_file_path)[0]) + '_regions.txt'
    autosomal_reference_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                                 'chr20', 'chr21', 'chr22'] # This should be generalized for non-human genomes
    input_fasta_file_path = os.path.abspath(input_fasta_file_path)
    
    if output_fasta_file_path is None:
        if only_autosomal:
            output_fasta_file_path = str(
                    os.path.splitext(input_fasta_file_path)[0]
                    ) + 'output_regions_only_autosomal.fa'
        else:
            output_fasta_file_path = str(os.path.splitext(input_fasta_file_path)[0]) + 'output_regions.fa'
    
    if regions is None and not only_autosomal:
        """
        If we are not filtering out any specific regions, just test for a fai file
        and create one if not found. Return the output_fasta_file_path.
        """
        if not os.path.isfile(input_fasta_file_path + '.fai'):
            # This always creates the index in the directory where the fasta file lives.
            # Don't create a copy of the reference if we're not filtering anything!
            faidx(input_fasta_file_path, '-o', output_fasta_file_path)
        return input_fasta_file_path

    # Otherwise, we need to filter out sequences we want and reindex the resulting fasta.
    if regions is None and only_autosomal:
        # Explicitly providing regions will always supercede only_autosomal
        if chrom_sizes == None:
            # We need chrom.sizes because the regions file must have start and end coordinates. See http://www.htslib.org/doc/samtools-faidx.html
            raise Exception("ERROR: create_fasta_file only_autosomal requires chrom_sizes.")
            return
        with open(regions_file, 'w') as rg_file:
            with open(chrom_sizes, 'r') as cs_file:
                for line in cs_file:
                    rec = line.split()
                    if rec[0] in autosomal_reference_names:
                        # Write a line to the regions file.
                        rg_file.write("{}:1-{}\n".format(rec[0], rec[1]))

    elif os.path.isfile(regions):
        regions_file = regions
        
    elif isinstance(regions, list):
        if only_autosomal:
            regions_list = []
            for region in regions:
                chrom = region.split(':')[0]
                if chrom in autosomal_reference_names:
                    regions_list.append(region)                
        else:
          regions_list = regions
            
        # Validate all entries in regions_list
        chrom_sizes_dict = {}
        with open(chrom_sizes, 'r') as cs_file:
            for line in cs_file:
                rec = line.split()
                chrom_sizes_dict[rec[0]] = rec[1]

        with open(regions_file, 'w') as file:
            for region in regions_list:
                # Quick and dirty validation -- just checks for presence of the proper delimiters
                # and end coordinate within the extent of the chromosome.
                s1 = region.split(':')
                s2 = region.split('-')
                if len(s1) == 2 and len(s2) == 2 and s1[0] in chrom_sizes_dict.keys() and s2[1] <= chrom_sizes_dict[s1[0]]:
                    file.write(str(region) + '\n')
                else:
                    if s1[0] in chrom_sizes_dict.keys() and s2[1] <= chrom_sizes_dict[s1[0]]:
                        sys.stderr.write('WARNING: Region {} coordinates not found in the reference sequence. Omitting.\n'.format(region))
                    else:
                        sys.stderr.write('WARNING: Region {} is not in the proper format (chrom:start-end). Omitting.\n'.format(region))

    faidx(input_fasta_file_path, '-r', regions_file, '-o', output_fasta_file_path)
    
    return output_fasta_file_path


def true_alignment_match(
        aligned_segment: AlignedSegment, needs_liftover: bool = False, contig: str = None, ref_start: int = None,
        ref_end: int = None, strand: str = None, tag_read: bool = False, only_output_match_label: bool = False, only_output_overlap: bool = False,
        true_reference_sequence_path: str = None, aligned_reference_sequence_path: str = None, liftover_converter = None,
        chain_file_path: str = None, sample: str = None, haplotype: int = None, reference_liftover_converters = None
        ) -> object:
    """

    Args:
        chain_file_path (object): 
        aligned_segment: 
        sample: 
        tag_read (object): 
    """
    overlap = 0
    match_label = 'non_match'
    
    if contig is None or ref_start is None or ref_end is None or strand is None:
        contig, ref_start, ref_end, strand = true_alignment_coordinates(
                aligned_segment, tag_read = tag_read, contig = contig, ref_start = ref_start, ref_end = ref_end,
                strand = strand
                )
    
    if needs_liftover:
        if liftover_converter is not None:
            liftover_converter = liftover_converter
    
        elif reference_liftover_converters is not None:
            if aligned_reference_sequence_path is not None:
                if aligned_reference_sequence_path in reference_liftover_converters:
                    if sample in reference_liftover_converters[aligned_reference_sequence_path]:
                        if str(haplotype) in reference_liftover_converters[aligned_reference_sequence_path][sample]:
                            liftover_converter = reference_liftover_converters[aligned_reference_sequence_path][sample][
                                str(haplotype)]
    
        elif chain_file_path is not None:
            if os.path.isfile(chain_file_path):
                liftover_converter = pyliftover.Liftover(chain_file_path)
    
        else:
            sys.stderr.write('Could not liftover\n')
            return
    
        contig = liftover_converter.convert_coordinate(contig, ref_start)[0][0]
        ref_start = liftover_converter.convert_coordinate(contig, ref_start)[0][1]
        ref_end = liftover_converter.convert_coordinate(contig, ref_end)[0][1]
        strand = liftover_converter.convert_coordinate(contig, ref_start)[0][2]
    
    if aligned_segment.reference_name == contig:
        match_label = 'ref_match'
        if int(ref_start) <= int(aligned_segment.reference_start) <= int(ref_end):
            match_label = 'mapping_match'
            if int(aligned_segment.reference_end) >= int(ref_end):
                overlap = int(ref_end) - int(aligned_segment.reference_start)
            else:
                overlap = int(aligned_segment.reference_end) - int(aligned_segment.reference_start)
        elif int(ref_start) <= int(aligned_segment.reference_end) <= int(ref_end):
            match_label = 'mapping_match'
            if int(ref_start) >= int(aligned_segment.reference_start):
                overlap = int(aligned_segment.reference_end) - int(ref_start)
            else:
                overlap = int(aligned_segment.reference_end) - int(aligned_segment.reference_start)
        elif (
                int(ref_start) - (int(ref_start) * 0.1)) <= aligned_segment.reference_start <= (
                int(ref_end) + (int(ref_end) * 0.1)):
            match_label = 'within_10percent'
    
    if tag_read:
        aligned_segment.set_tag(
                tag = 'ov', value = str(overlap / (int(ref_end) - int(ref_start))), value_type = 'Z', replace = True
                )
        aligned_segment.set_tag(
                tag = 'OV', value = str(overlap / (aligned_segment.query_alignment_length)), value_type = 'Z',
                replace = True
                )
        aligned_segment.set_tag(tag = 'ml', value = str(match_label), value_type = 'Z', replace = True)
    
    if only_output_match_label:
        return match_label
    
    elif only_output_overlap:
        return overlap
    
    else:
        return match_label, overlap


def true_alignment_coordinates(
        aligned_segment: AlignedSegment, tag_read: bool = False, contig: str = None, ref_start: int = None,
        ref_end: int = None, strand: str = None
        ) -> object:
    """

    Args:
        aligned_segment: 
        tag_read: 
        contig: 
        ref_start: 
        ref_end: 
        strand: 

    Returns:

    """
    if '!' in str(aligned_segment.query_name):
        _read_name, _contig, _ref_start, _ref_end, _strand = aligned_segment.query_name.split('!')
        if contig is None:
            contig = _contig
        if ref_start is None:
            ref_start = _ref_start
        if ref_end is None:
            ref_end = _ref_end
        if strand is None:
            strand = _strand
    
    if tag_read:
        aligned_segment.set_tag(tag = 'st', value = str(ref_start), value_type = 'Z', replace = True)
        aligned_segment.set_tag(tag = 'lo', value = str(contig), value_type = 'Z', replace = True)
        aligned_segment.set_tag(tag = 'en', value = str(ref_end), value_type = 'Z', replace = True)
        aligned_segment.set_tag(tag = 'sd', value = str(strand), value_type = 'Z', replace = True)
    
    return contig, ref_start, ref_end, strand


def alignment_type(aligned_segment: AlignedSegment, tag_read: bool = False) -> str:
    """

    Args:
        aligned_segment:
        tag_read:

    Returns:

    """
    alignment_label: str = 'None'
    if aligned_segment.is_unmapped:
        alignment_label = 'unmapped'
    elif aligned_segment.is_secondary:
        alignment_label = 'secondary'
    elif aligned_segment.is_supplementary:
        alignment_label = 'supplementary'
    else:
        alignment_label = 'mapped'
    if tag_read:
        aligned_segment.set_tag(tag = 'al', value = str(alignment_label), value_type = 'Z', replace = True)
    
    return alignment_label

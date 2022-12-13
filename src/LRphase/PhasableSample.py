# coding=utf-8
import time
from typing import Any, Iterator
import sys, os
import pysam
from LRphase import SimulatePhasedData
import pyliftover


class PhasableSample:
    """

    """
    #RG_ID_dict: object 
    
    def __init__(
            self,
            sample,
            vcf_file_path,
            ignore_phase_sets,
            alignment_file_paths,
            RG_ID_dict,
            bam_header = None,
            reference_sequence_names = None,
            reference_sequence_paths = None,
            simulate_haplotypes = False,
            output_directory_path = None,
            only_autosomal = False,
            regions = None,
            ploidy = 2,
            quiet = False,
            silent = False
    ):
        self.sample = sample
        
        self.vcf_file = pysam.VariantFile(vcf_file_path)
        self.vcf_file_path = vcf_file_path
        self.ignore_phase_sets = ignore_phase_sets
        
        self.alignment_file_paths = []
        self.RG_ID_dict = RG_ID_dict
        self.use_RG_tag = {}
        
        #sys.stderr.write("%s\n" %(alignment_file_paths))        
        for alignment_file_path in alignment_file_paths:
            self.alignment_file_paths.append(str(alignment_file_path))
            self.use_RG_tag[str(alignment_file_path)] = alignment_file_paths[str(alignment_file_path)]

        if bam_header:
            self.bam_header = bam_header
        else:
            # Take the header from the first alignment file.
            self.bam_header = pysam.AlignmentFile(self.alignment_file_paths[0], 'rb').header

        self.reference_sequence_names = reference_sequence_names
        self.reference_sequence_paths = reference_sequence_paths
        self.reference_sequences_in_VCF = list(self.vcf_file.index)
        self.only_autosomal = only_autosomal
        self.regions = regions
        self.ploidy = ploidy
        self.quiet = quiet
        self.silent = silent
        
        if output_directory_path is None:
            self.output_directory = str(self.sample)
        else:
            self.output_directory = output_directory_path
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)
        
        if simulate_haplotypes:
            
            self.reference_sequence_path = str(os.path.abspath(self.output_directory)) + '/reference_sequence.fa'
            reference_sequence_list = []
            for i, reference_sequence in enumerate(self.reference_sequence_paths):
                reference_sequence_list.append(SimulatePhasedData.create_fasta_file(input_fasta_file_path = reference_sequence, output_fasta_file_path = str(os.path.abspath(self.output_directory)) + '/reference_sequence'+str(i)+'.fa', regions = self.regions, only_autosomal = self.only_autosomal))
            
            with open(self.reference_sequence_path, 'a') as outfile:
                for file in reference_sequence_list:
                    with open(file, 'r') as infile:
                        outfile.write(infile.read())
                        
            SimulatePhasedData.index_fasta_file(self.reference_sequence_path)
            
            self.haplotype_reference_sequences = {}
            self.haplotype_chain_files = {}
            self.haplotype_liftover_converters = {}
            for j in range(1,self.ploidy+1):
                self.haplotype_reference_sequences[j], self.haplotype_chain_files[j] = SimulatePhasedData.generate_haplotype_specific_fasta(j, self.sample, self.reference_sequence_path, self.vcf_file_path, output_reference_sequence_path = str(os.path.abspath(self.output_directory)) + '/haplotype'+str(j)+'_reference_sequence.fa', chain_file_path = str(os.path.abspath(self.output_directory)) + '/haplotype'+str(j)+'_reference_sequence.chain')
                self.haplotype_liftover_converters[j] = pyliftover.LiftOver(self.haplotype_chain_files[j])

    def _pysam_bam_files_initialize(self):
        bam_files = []
        for alignment_file_path in self.alignment_file_paths:
            bam_files.append(iter(pysam.AlignmentFile(alignment_file_path, 'rb')))
        return bam_files

    def __iter__(self):
        #sys.stderr.write('Alignments for sample %s\n' % str(self.sample))
        self.bam_files = self._pysam_bam_files_initialize()
        # self._initialize_alignment_counter()
        self.alignment_files_processed_count = 0
        self.alignments_processed_count = 0
        self.alignment_files_read_counts = []
        bam_file: Iterator[Any]
        for bam_file in self.bam_files:
            self.alignment_files_read_counts.append(bam_file.mapped + bam_file.unmapped)
        self.total_alignments = sum(self.alignment_files_read_counts)
        self.alignment_file_pysam = self.bam_files[self.alignment_files_processed_count]
        self.start_process_time = time.time()
        self.unique_read_names = set()
        return self

    def __next__(self):
        read = next(self.alignment_file_pysam)
        if read.query_name not in self.unique_read_names:
            self.unique_read_names.add(read.query_name)
        #read = self._evaluate_alignment(read)
        RG_info: object = self._get_RG_info_for_read(read, self.alignment_file_pysam)
        #if str(RG_info[2])[0:3] == 'SIM':
        #read.set_tag(tag='oa', value=str(RG_info[2][3]), value_type='Z', replace=True)
        read.set_tag(tag = 'RG', value = str(RG_info[1]), value_type = 'Z', replace = True)
        self.alignment_files_read_counts[self.alignment_files_processed_count] -= 1
        if self.alignment_files_read_counts[self.alignment_files_processed_count] == 0:
            self.alignment_files_processed_count += 1
            if self.alignment_files_processed_count == len(self.alignment_file_paths):
                raise StopIteration()
            self.alignment_file_pysam = self.bam_files[self.alignment_files_processed_count]
        self.alignments_processed_count += 1
        if self.alignments_processed_count % 2000 == 0 and not (self.quiet or self.silent):
            print(
                    'Processing %s reads per second' % str(
                            round(self.alignments_processed_count / (time.time() - self.start_process_time), 2)
                            )
                    )
            print(
                    'Processed %s reads (%s percent complete)' % (str(self.alignments_processed_count), str(
                            round(100 * (self.alignments_processed_count / self.total_alignments), 2)
                            ))
                    )
            print(
                    'Processing sample %s will finish in %s seconds' % (str(self.sample), str(
                            round(
                                    (self.total_alignments - self.alignments_processed_count) / (
                                            self.alignments_processed_count / (time.time() - self.start_process_time)),
                                    1
                                    )
                            ))
                    )
        return read

    def __repr__(self):
        return f'Sample Name: {self.sample}\nVCF: {self.vcf_file_path}\nReference sequence path: {str(self.reference_sequence_paths)}\nTotal sequences in reference files: {len(self.reference_sequence_names)}\nTotal Reference sequences in VCF: {len(self.reference_sequences_in_VCF)}\nAlignment Files: {str(self.alignment_file_paths)}\nTotal alignment files processed: {self.alignment_files_processed_count}\nTotal alignments: {self.total_alignments}\nTotal alignments processed: {self.alignments_processed_count}\nTotal unique reads observed: {len(self.unique_read_names)}'

    def _get_RG_info_for_read(self, read: pysam.AlignedSegment, alignment_file: object) -> object:
        """

        Args:
            read: 
            alignment_file: 

        Returns:

        """
        alignment_file_path = ''
        if isinstance(alignment_file, pysam.AlignmentFile):
            if str(alignment_file.filename.decode()) in self.alignment_file_paths:
                alignment_file_path = str(alignment_file.filename.decode())
            else:
                return
        elif str(alignment_file) in self.alignment_file_paths:
            alignment_file_path = str(alignment_file)
        else:
            return
        if alignment_file_path in self.use_RG_tag.keys():
            use_RG_tag = self.use_RG_tag[alignment_file_path]
            if use_RG_tag:
                if str(read.get_tag('RG')) in self.RG_ID_dict[alignment_file_path]:
                    if self.RG_ID_dict[alignment_file_path][str(read.get_tag('RG'))]['SM'] == self.sample:
                        ID = str(read.get_tag('RG'))
                        outputID = self.RG_ID_dict[alignment_file_path][str(read.get_tag('RG'))]['outputID']
                        sample_description = self.RG_ID_dict[alignment_file_path][str(read.get_tag('RG'))]['DS']
                        return ID, outputID, sample_description, use_RG_tag
            else:
                ID = list(self.RG_ID_dict[alignment_file_path].keys())[0]
                outputID = self.RG_ID_dict[alignment_file_path][ID]['outputID']
                sample_description = self.RG_ID_dict[alignment_file_path][ID]['DS']
                return ID, outputID, sample_description, use_RG_tag

            
    def simulate_reads(self, path_to_pbsim = 'pbsim', depth = 1,
        simulation_mode = 'pbsim2/data/R103.model',
        difference_ratio = '23:31:46', length_mean = 20000,
        length_max = 1000000, length_min = 100,
        length_sd = 15000, accuracy_min = 0.01,
        accuracy_max = 1.00, accuracy_mean = 0.80,
        prefix = None, id_prefix = 'S', output_directory = None, sample = None,
        haplotypes = [1,2]):
        
        self.haplotype_simulated_reads = {}
        for haplotype in haplotypes:
            self.haplotype_simulated_reads[haplotype] = SimulatePhasedData.simulate_reads_pbsim2(self.haplotype_reference_sequences[haplotype],
                                                                                                 path_to_pbsim = path_to_pbsim,
                                                                                                 depth = depth,
                                                                                                 simulation_mode = simulation_mode,
                                                                                                 difference_ratio = difference_ratio,
                                                                                                 length_mean = length_mean,
                                                                                                 length_max = length_max,
                                                                                                 length_min = length_min,
                                                                                                 length_sd = length_sd,
                                                                                                 accuracy_min = accuracy_min,
                                                                                                 accuracy_max = accuracy_max,
                                                                                                 accuracy_mean = accuracy_mean,
                                                                                                 prefix = prefix,
                                                                                                 id_prefix = id_prefix,
                                                                                                 output_directory = self.output_directory +
                                                                                                 '/simulated',
                                                                                                 sample = self.sample,
                                                                                                 haplotype = haplotype
            )
            
            
    
        
                       
        
    # @classmethod
    # def PhasableSample(
    #         cls, sample, vcf_file_path, ignore_phase_sets, param, RG_ID_dict: object, reference_sequence_names,
    #         reference_sequence_paths
    #         ):
    #     """
    #
    #     Args:
    #         sample:
    #         vcf_file_path:
    #         ignore_phase_sets:
    #         reference_sequence_names:
    #         reference_sequence_paths:
    #         RG_ID_dict (object):
    #
    #     Returns:
    #         object:
    #     """
    #     pass

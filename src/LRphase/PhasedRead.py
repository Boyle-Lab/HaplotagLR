import math
import sys, os
from collections import defaultdict
from pysam import VariantFile, AlignmentFile
#import numpy as np
import pyliftover
from LRphase import PhasableSample, PhaseSet


def true_alignment_match(
    aligned_segment: object,
    needs_liftover: object = False,
    contig: object = None,
    ref_start: object = None,
    ref_end: object = None,
    strand: object = None,
    tag_read: object = False,
    only_output_match_label: object = False,
    only_output_overlap: object = False,
    true_reference_sequence_path: object = None,
    aligned_reference_sequence_path: object = None,
    liftover_converter: object = None,
    chain_file_path: object = None,
    sample: object = None,
    haplotype: object = None,
    reference_liftover_converters: object = None
) -> object:
    """
    Args:
    aligned_segment:
    needs_liftover:
    contig:
    ref_start (object):
    ref_end:
    strand:
    tag_read:
    only_output_match_label:
    only_output_overlap:
    true_reference_sequence_path:
    aligned_reference_sequence_path (object):
    liftover_converter:
    chain_file_path:
    sample:
    haplotype:
    reference_liftover_converters:
    
    Returns:
    
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
            sys.stderr.write('could not liftover\n')
            return

        ref_start_lifted = liftover_converter.convert_coordinate(contig, ref_start)
        ref_end_lifted = liftover_converter.convert_coordinate(contig, ref_end)
        
        contig = ref_start_lifted[0][0]
        ref_start = ref_start_lifted[0][1]
        ref_end = ref_end_lifted[0][1]
        strand = ref_start_lifted[0][2]

    """
    sys.stderr.write("True: %s:%s-%s, %s; Observed: %s:%s-%s\n" % (contig, ref_start, ref_end, strand,
                                                                   aligned_segment.reference_name,
                                                                   aligned_segment.reference_start,
                                                                   aligned_segment.reference_end))
    """
    if aligned_segment.reference_name == contig:
        match_label = 'ref_match'
        if int(aligned_segment.reference_start) >= int(ref_start) and int(aligned_segment.reference_start) <= int(
                ref_end
                ):
            match_label = 'mapping_match'
            if int(aligned_segment.reference_end) >= int(ref_end):
                overlap = int(ref_end) - int(aligned_segment.reference_start)
            else:
                overlap = int(aligned_segment.reference_end) - int(aligned_segment.reference_start)
        elif int(aligned_segment.reference_end) >= int(ref_start) and int(aligned_segment.reference_end) <= int(
                ref_end
                ):
            match_label = 'mapping_match'
            if int(ref_start) >= int(aligned_segment.reference_start):
                overlap = int(aligned_segment.reference_end) - int(ref_start)
            else:
                overlap = int(aligned_segment.reference_end) - int(aligned_segment.reference_start)
        elif aligned_segment.reference_start >= (
                int(ref_start) - (int(ref_start) * 0.1)) and aligned_segment.reference_start <= (
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
        aligned_segment, tag_read = False, contig = None, ref_start = None, ref_end = None, strand = None
        ):
    if '!' in str(aligned_segment.query_name):
        _read_name, _contig, _ref_start, _ref_end, _strand = aligned_segment.query_name.split('!')
        if contig is None:
            contig = _contig.split(':')[0]
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


def true_read_origin(
        aligned_segment, tag_read = False, sample = None, haplotype = None, contig = None, ref_start = None,
        ref_end = None, strand = None
        ):
    if '!' in str(aligned_segment.query_name):
        _ref_start: object
        _read_name, _contig, _ref_start, _ref_end, _strand = aligned_segment.query_name.split('!')
        if contig is None:
            contig = _contig
        if ref_start is None:
            ref_start = _ref_start
        if ref_end is None:
            ref_end = _ref_end
        if strand is None:
            strand = _strand
        
        if '__' in _read_name:
            sample, haplotype, name_contig_numeric = _read_name.split('__')
    
    if tag_read:        
        aligned_segment.set_tag(tag = 'sm', value = str(sample), value_type = 'Z', replace = True)
        aligned_segment.set_tag(tag = 'ha', value = str(haplotype), value_type = 'Z', replace = True)
    
    return sample, haplotype


def alignment_type(aligned_segment, tag_read = False):
    alignment_type = 'None'
    if aligned_segment.is_unmapped:
        alignment_type = 'unmapped'
    elif aligned_segment.is_secondary:
        alignment_type = 'secondary'
    elif aligned_segment.is_supplementary:
        alignment_type = 'supplementary'
    else:
        alignment_type = 'mapped'
    if tag_read:
        aligned_segment.set_tag(tag = 'al', value = str(alignment_type), value_type = 'Z', replace = True)
    
    return alignment_type


class PhasedRead:
    """

    """
    
    def __init__(
            self, aligned_segment, phasable_sample = None, vcf_file = None, sample = None,
            haplotype = None, ignore_phase_sets = False, error_model = 0,
            error_rate_threshold = 0.01,
            prior_probabilities = None, bam_file_header = None, bam_file_template = None,
            output_file_path = None, liftover_converters = None,
            multinomial_correction = True, auto_phase = True, evaluate_alignment = True,
            evaluate_true_alignment = False, aligned_reference_sequence_path = None,
            powlaw_alpha = 4.5, powlaw_xmin = 2.0
    ):
        self.output_file_path = None
        self.liftover_converters = None
        self.aligned_segment = aligned_segment
        self.aligned_reference_sequence_path = None
        self.needs_liftover = False
        self.true_reference_sequence_path = None
        
        if isinstance(phasable_sample, PhasableSample.PhasableSample):
            self.vcf_file=phasable_sample.vcf_file
            self.sample=phasable_sample.sample
            if phasable_sample.simulate_haplotypes:
                if output_file_path is None:
                    self.output_file_path = phasable_sample.output_directory
                self.liftover_converters = phasable_sample.haplotype_liftover_converters
                self.aligned_reference_sequence_path = phasable_sample.reference_sequence_path
                self.needs_liftover = True
                self.evaluate_true_alignment = True
        
        else:
            self.vcf_file = vcf_file
            self.sample = sample
        
        self.haplotype = haplotype
        self.ignore_phase_sets = ignore_phase_sets

        # Attributes related to error-rate thresholding
        self.error_model = error_model
        self.multinomial_correction = multinomial_correction
        self.error_rate_threshold = error_rate_threshold
        self.prior_probabilities = prior_probabilities
        self.powlaw_alpha = powlaw_alpha
        self.powlaw_xmin = powlaw_xmin

        # self._get_alignment_label()
        if self.liftover_converters is None:
            self.liftover_converters = liftover_converters
        self.auto_phase = auto_phase
        self._PhaseSet_max = None
        self.PhaseSets = []
        self.evaluate_alignment = evaluate_alignment
        self.evaluate_true_alignment = evaluate_true_alignment
        
        if output_file_path:
            self.output_file_path = output_file_path
        else:
            self.output_file_path = '%s_phase_tagged.bam' % self.sample
            
        if bam_file_header:
            self.bam_file_header = bam_file_header
        else:
            self.bam_file_header = None

        if bam_file_template:
            self.bam_file_template = AlignmentFile(bam_file_template, 'rb')
        else:
            self.bam_file_template = None
            
        if self.evaluate_alignment:
            alignment_type(self.aligned_segment, tag_read = True)
            #self._evaluate_alignment()
        if self.evaluate_true_alignment:
            true_alignment_match(
                self.aligned_segment, 
                needs_liftover = self.needs_liftover, 
                contig = None,
                ref_start = None,
                ref_end = None,
                strand = None,
                tag_read = True,
                only_output_match_label = False,
                only_output_overlap = False,
                true_reference_sequence_path = self.true_reference_sequence_path,
                aligned_reference_sequence_path = self.aligned_reference_sequence_path,
                liftover_converter = self.liftover_converters,
                chain_file_path = None,
                sample = self.sample,
                haplotype = self.haplotype,
                reference_liftover_converters = None
            )
            true_read_origin(
                self.aligned_segment,
                tag_read = True,
                sample = self.sample,
                haplotype = self.haplotype,
                contig = None,
                ref_start = None,
                ref_end = None,
                strand = None
            )
        
        if self.auto_phase:
            self.phase_read(
                    error_model = self.error_model, error_rate_threshold = self.error_rate_threshold,
                    multinomial_correction = self.multinomial_correction,
                    )
    
    def __repr__(self):
        return f'Read Name: {self.query_name}, Is phased: {self.is_phased}, Favors single phase: {self.one_phase_is_favored}\nRead Length: {self.aligned_segment.query_length}\nAlignment Type: {self.alignment_type}\nPrimary Alignment Length: {self.aligned_segment.query_alignment_length}'
    
    def __iter__(self):
        self.PhaseSets = self._find_PhaseSets(
                self.error_model, self.error_rate_threshold, self.prior_probabilities,
                liftover_converters = self.liftover_converters
                )
        self.PhaseSets_processed_count = 0
        # self.alignment_files_read_counts = [bam_file.mapped+bam_file.unmapped for bam_file in self.bam_files]
        return self
    
    def __next__(self):
        if self.PhaseSets_processed_count < len(self.PhaseSets):
            PhaseSet = self.PhaseSets[self.PhaseSets_processed_count]
            self.PhaseSets_processed_count += 1
            return PhaseSet
        else:
            raise StopIteration()
    
    def _find_PhaseSets(self, error_model, error_rate_threshold, prior_probabilities, liftover_converters = None):
        
        if liftover_converters is None:
            liftover_converters = self.liftover_converters
        
        al_type = alignment_type(self.aligned_segment, tag_read=True)
        
        if self.aligned_segment.is_unmapped:
            # Return an empty phaseSet
            return []
        
        if self.aligned_segment.get_tag('al') == 'aligned_to_contig_not_in_vcf' or self.aligned_segment.is_unmapped:
            return 'aligned_to_contig_not_in_vcf'
        
        variants_phase_sets = self._fetch_phased_variants(
                self.aligned_segment, self.vcf_file, self.sample, self.ignore_phase_sets
                )
        
        PhaseSets = []
        if len(variants_phase_sets) > 0:
            if max([len(variants_phase_sets[str(phase_set)]) for phase_set in variants_phase_sets]) > 0:
                for phase_set in variants_phase_sets:
                    PhaseSets.append(
                        PhaseSet.PhaseSet(
                            phase_set,
                            self.aligned_segment,
                            variants_phase_sets[str(phase_set)],
                            str(self.sample),
                            error_model = error_model,
                            error_rate_threshold = error_rate_threshold,
                            prior_probabilities = prior_probabilities,
                            liftover_converters = liftover_converters,
                            powlaw_alpha = self.powlaw_alpha,
                            powlaw_xmin = self.powlaw_xmin
                        )
                    )
        
        return PhaseSets
    
    def phase_read(
        self, 
        error_model = None,
        error_rate_threshold = None, 
        prior_probabilities = None,
        multinomial_correction = None
    ):
        '''#self.PhaseSet_max, self.phase, self.phase_set_name, self.log_likelihood_ratio, self.max_phase,
        self.ploidy_phase_set = self._select_best_phasing(self._find_PhaseSets()) '''
        if error_model is None:
            error_model = self.error_model
        if error_rate_threshold is None:
            error_rate_threshold = self.error_rate_threshold
        if multinomial_correction is None:
            multinomial_correction = True
        
        #if evaluate_alignment:
            #self._evaluate_alignment()
        #if evaluate_true_alignment:
        #    if self.aligned_segment.has_tag('oa'):
        #        match_label = self._evaluate_true_alignment
        
        if not self.aligned_segment.has_tag('al'):
            alignment_type(self.aligned_segment,tag_read=True)
        
        #if not self.aligned_segment.has_tag('ml'):
        #    if self.aligned_segment.has_tag('oa'):
        #        match_label = self._evaluate_true_alignment
        
        #_PhaseSet_max, phase, phase_set_name, log_likelihood_ratio, max_phase = self._select_best_phasing(
        _PhaseSet_max = self._select_best_phasing(
            self._find_PhaseSets(
                error_model,
                error_rate_threshold,
                prior_probabilities,
                liftover_converters = self.liftover_converters
            ),
            error_model,
            error_rate_threshold,
            prior_probabilities,
            multinomial_correction
        )
        
        self._PhaseSet_max = _PhaseSet_max
        
        self.aligned_segment.set_tag(tag = 'PS', value = str(self.phase_set_name), value_type='Z', replace=True)
        self.aligned_segment.set_tag(tag = 'HP', value = str(self.phase), value_type='Z', replace=True)
        self.aligned_segment.set_tag(tag = 'PC', value = str(self.log_likelihood_ratio), value_type='Z', replace=True)
        self.aligned_segment.set_tag(tag = 'py', value = str(self.ploidy_phase_set), value_type='Z', replace=True)
        self.aligned_segment.set_tag(tag = 'HS', value = str(self._PhaseSet_max.total_hets), value_type='Z', replace=True)
        
        return self
    
    def _select_best_phasing(
            self, PhaseSets, error_model, error_rate_threshold, prior_probabilities = None,
            multinomial_correction = True
            ):
        '''self.PhaseSet_max, self.phase, self.phase_set_name, self.log_likelihood_ratio, self.max_phase =
        self._select_best_phasing(self.PhaseSets) '''
        self.PhaseSets = PhaseSets
        phasable = {}
        if isinstance(PhaseSets, str):
            return PhaseSet.PhaseSet(None, None, None, None)
        
        if len(PhaseSets) > 0:
            for _PhaseSet in PhaseSets:
                _PhaseSet.solve_phase(
                    error_model,
                    error_rate_threshold,
                    prior_probabilities,
                    multinomial_correction
                )
                if isinstance(_PhaseSet.max_log_likelihood_ratio, float):
                    phasable[_PhaseSet.max_log_likelihood_ratio] = _PhaseSet
        
        if len(phasable) > 0:
            PhaseSet_max = phasable[max(phasable.keys())]
        else:
            PhaseSet_max = PhaseSet.PhaseSet(None, self.aligned_segment, [], str(self.sample))
        
        return PhaseSet_max#, phase, phase_set_name, log_likelihood_ratio, max_phase
    
    # @staticmethod
    def _fetch_phased_variants(self, aligned_segment, vcf_file, sample, ignore_phase_sets):
        '''self.variants_phase_sets = self._fetch_phased_variants(self.aligned_segment, self.vcf_file, self.sample,
        self.ignore_phase_sets) '''
        variants_phase_sets = defaultdict(list)

        vcf_in = VariantFile(vcf_file)

        # VariantFile.fetch() will cause program termination if it encounters a contig
        # that is not found in the given VCF file. Therefore we need to encapsulate this
        # in a try block to prevent premature termination. Aligned segments that map to
        # such contigs will be skipped.
        try:
            vcf_recs = vcf_in.fetch(aligned_segment.reference_name, aligned_segment.reference_start, aligned_segment.reference_end)
            # Note that coordinates supplied in this way are assumed zero-based by pysam (see source code, https://github.com/pysam-developers/pysam/blob/master/pysam/libcbcf.pyx, line 4366. We should verify these are zero-based. If not, use pysam chr:start-end string -- these are treated as one-based. This may not be a huge deal, since only variants coinciding with the start/end of the fragment will be affected.
            
            for vcf_record in vcf_recs:
                if vcf_record.samples[str(sample)].phased and max(vcf_record.samples[str(sample)].allele_indices) > min(
                        vcf_record.samples[str(sample)].allele_indices
                        ):
                    if ignore_phase_sets:
                        variants_phase_sets['ignore_phase_sets'].append(vcf_record)
                    else:
                        variants_phase_sets[str(vcf_record.samples[str(sample)]['PS'])].append(vcf_record)
        except:
            # Do some error handling...
            pass
        
        vcf_in.close()
        return variants_phase_sets

    def write_to_bam(self, output_file_path_bam=None, output_bam_pysam=None, header=None, template=None):
        """
        Args:
            output_file_path_bam: Path to output bam file.
            output_bam_pysam: Pysam AlignmentFile object invoked with 'wb' flags.
            header: Pysam AlignmentHeader object or dictionary.
            template: Pysam AlignmentFile object invoked with 'rb' flags from which bam header will be extracted.
        Returns: Nothing
        """

        if header and template:
            sys.stderr.write('WARNING: Both header and template given for bam output. Template will be used by default.')        
        _bam_file_header = None
        if header:
            _bam_file_header = header
        elif self.bam_file_header:
            _bam_file_header = self.bam_file_header

        _bam_file_template = None
        if template:
            _bam_file_template = template
        elif self.bam_file_template:
            _bam_file_template = self.bam_file_template
            
        if not _bam_file_header and not _bam_file_template and not output_bam_pysam:
            raise Exception('Bam output requires header, template, or pysam AlignmentFile object (output_bam_pysam). None given.')
            
        if output_bam_pysam:
            output_bam = output_bam_pysam
        else:
            if output_file_path_bam:
                _output_file_path = output_file_path_bam
            else:
                _output_file_path = self.output_file_path
            if _bam_file_template:
                output_bam = AlignmentFile(_output_file_path, 'wb', template = _bam_file_template)
            else:
                output_bam = AlignmentFile(_output_file_path, 'wb', header = _bam_file_header)

        # self.aligned_segment.set_tag(tag = 'RG', value = str(RG_info[1]), value_type='Z', replace=True)
        output_bam.write(self.aligned_segment)

        return

    
    @property
    def is_phased_correctly(self):
        '''
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            if self.is_phased:
                return False
            else:
                return
        else:
            if self.is_phased:
                if self.aligned_segment.has_tag('ha'):
                    return str(self.phase) == str(self.aligned_segment.get_tag('ha'))
    
    @property
    def read_name(self):
        '''
        '''
        return str(self.aligned_segment.query_name)
    
    @property
    def query_name(self):
        '''
        '''
        return str(self.aligned_segment.query_name)
    
    @property
    def alignment_type(self):
        '''
        '''
        if self.aligned_segment.has_tag('al'):
            return str(self.aligned_segment.get_tag('al'))
    
    @property
    def fraction_of_read_overlapping_true_alignment(self):
        '''
        '''
        if self.aligned_segment.has_tag('ov'):
            return str(self.aligned_segment.get_tag('ov'))
    
    @property
    def alignment_is_unmapped(self):
        '''
        '''
        return self.aligned_segment.is_unmapped
    
    @property
    def alignment_is_supplementary(self):
        '''
        '''
        return self.aligned_segment.is_supplementary
    
    @property
    def alignment_is_secondary(self):
        '''
        '''
        return self.aligned_segment.is_secondary
    
    @property
    def is_aligned_to_contig_not_in_vcf(self):
        '''
        '''
        return self.aligned_segment.reference_name not in list(self.vcf_file.index)
    
    @property
    def alignment_is_reverse(self):
        '''
        true if aligned_segment is mapped to reverse strand
        '''
        return self.aligned_segment.is_reverse
    
    @property
    def alignment_is_mapped(self):
        '''
        '''
        return self.aligned_segment.is_unmapped == False
    
    @property
    def alignment_is_primary(self):
        '''
        '''
        if self.aligned_segment.has_tag('tp'):
            return str(self.aligned_segment.get_tag('tp')) == 'P'
        # return self.aligned_segment.is_unmapped == False
    
    @property
    def validation_using_true_alignment_label(self):
        '''
        '''
        if self.aligned_segment.has_tag('ml'):
            return str(self.aligned_segment.get_tag('ml'))
    
    @property
    def matches_true_alignment(self):
        '''
        '''
        if self.aligned_segment.has_tag('ml'):
            return str(self.aligned_segment.get_tag('ml')) == 'mapping_match'
    
    @property
    def mapped_within_10percent_of_true_alignment(self):
        '''
        '''
        if self.aligned_segment.has_tag('ml'):
            return str(self.aligned_segment.get_tag('ml')) == 'within_10percent'
    
    @property
    def does_not_overlap_true_alignment_but_does_align_to_true_contig(self):
        '''
        '''
        if self.aligned_segment.has_tag('ml'):
            return str(self.aligned_segment.get_tag('ml')) == 'ref_match'
    
    @property
    def does_not_match_true_alignment(self):
        '''
        '''
        if self.aligned_segment.has_tag('ml'):
            return str(self.aligned_segment.get_tag('ml')) == 'non_match'
    
    @property
    def alignment_length_fraction_of_total_read_length(self):
        '''
        aligned_segment.query_alignment_length/aligned_segment.query_length
        returns 0 if unmapped
        '''
        # self.aligned_segment.set_tag(tag = 'xx', value = str(self.aligned_segment.query_alignment_length/self.aligned_segment.query_length),
        # value_type='Z', replace=True)
        if self.aligned_segment.has_tag('al'):
            if str(self.aligned_segment.get_tag('al')) == 'unmapped_reads':
                return 0
            return self.aligned_segment.query_alignment_length / self.aligned_segment.query_length
    
    @property
    def fraction_of_alignment_length_overlapping_true_alignment(self):
        '''
        '''
        if self.aligned_segment.has_tag('OV'):
            return str(self.aligned_segment.get_tag('OV'))
    
    @property
    def favored_phase_is_correct(self):
        '''
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            if self.one_phase_is_favored:
                if self.aligned_segment.has_tag('ha'):
                    return str(self.max_phase) == str(self.aligned_segment.get_tag('ha'))
    
    @property
    def is_phased_incorrectly(self):
        '''
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            if self.is_phased:
                if self.aligned_segment.has_tag('ha'):
                    return str(self.phase) != str(self.aligned_segment.get_tag('ha'))
    
    @property
    def favored_phase_is_incorrect(self):
        '''
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            if self.one_phase_is_favored:
                if self.aligned_segment.has_tag('ha'):
                    return str(self.max_phase) != str(self.aligned_segment.get_tag('ha'))
    
    @property
    def is_Nonphasable(self):
        '''
        Returns True if this aligned_segment was assigned Nonphasable for not having any heterozygous positions
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return True
        else:
            return self.phase == 'Nonphasable'
    
    @property
    def is_Unphased(self):
        '''
        Returns True if this aligned_segment was left unassigned (Unphased) because it did not have a log likelihood ratio (LR)
        above the threshold chosen for this experiment.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            return self.phase == 'Unphased'
    
    @property
    def is_phased(self):
        '''
        Returns True if this aligned_segment was assigned to a haplotype based on its log likelihood ratio (LR) being above the
        threshold chosen for this experiment.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            if self.phase == 'Unphased' or self.phase == 'Nonphasable':
                return False
            else:
                return True
        return False
    
    @property
    def one_phase_is_favored(self):
        '''
        Returns True if this aligned_segment was assigned to a haplotype based on its log likelihood ratio (LR) being above the
        threshold chosen for this experiment.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return False
        else:
            if int(self.max_phase) > 0:
                one_phase_is_favored = True
            elif int(self.max_phase) == 0:
                one_phase_is_favored = False
            else:
                return
            return one_phase_is_favored
    
    def is_assigned_to_haplotype_i(self, haplotype):
        
        '''Returns True if this aligned_segment was was assigned to haplotype i. (ie: the index of the haplotype in the allele_indices A|T )
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        else:
            return str(self.phase) == str(haplotype)
    
    @property
    def overlaps_multiple_phase_sets(self):
        '''
        This aligned_segment overlapped heterozygous positions that belonged to at least two different phase sets (independently
        phased blocks of sequence are only phased within themselves and one must be chosen)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max == 'Nonphasable':
            return False
        if len(self.PhaseSets) > 1:
            overlaps_multiple_phase_sets = True
        else:
            overlaps_multiple_phase_sets = False
        return overlaps_multiple_phase_sets
    
    @property
    def PhaseSet_max(self):
        '''
        '''
        if self._PhaseSet_max is None:
            return
        elif self._PhaseSet_max == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self._PhaseSet_max
    
    @property
    def phase(self):
        '''
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.phase
    
    @property
    def max_phase(self):
        '''
        The number of distinct haplotypes at this position. Our model assumes that ploidy represents the ground truth
        molecular copy number of the population from which these reads were samples.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.max_phase
    
    @property
    def phase_set_name(self):
        '''
        The number of distinct haplotypes at this position. Our model assumes that ploidy represents the ground truth
        molecular copy number of the population from which these reads were samples.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.phase_set
    
    @property
    def phasing_error_rate(self):
        if self.phase == 'Nonphasable':
            return 0.50
        if self.log_likelihood_ratio > 0:
            p = self.powlaw_modified(self.log_likelihood_ratio,
                                     a = self.powlaw_alpha,
                                     xmin = self.powlaw_xmin)
            pl = 1 - (1 / float(self.ploidy_phase_set))
            if p <= pl:
                return p
            else:
                return pl
    
    def powlaw_modified(self, x, a = 4.5, xmin = 2):
        return ((a - 1) / xmin) * math.pow(x / xmin, -1 * a)
    
    @property
    def log_likelihood_ratio(self):
        '''
        The number of distinct haplotypes at this position. Our model assumes that ploidy represents the ground truth
        molecular copy number of the population from which these reads were samples.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.max_log_likelihood_ratio
    
    @property
    def ploidy_phase_set(self):
        '''
        The number of distinct haplotypes at this position. Our model assumes that ploidy represents the ground truth
        molecular copy number of the population from which these reads were samples.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.ploidy_phase_set
    
    @property
    def error_model_used(self):
        '''
        self.aligned_segment.set_tag(tag = 'em', value = str(self.error_model), value_type='Z', replace=True)

        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.error_model_used
    
    @property
    def multinomial_correction_used(self):
        '''
        self.aligned_segment.set_tag(tag = 'em', value = str(self.error_model), value_type='Z', replace=True)

        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return False
        else:
            return self.PhaseSet_max.multinomial_correction
    
    @property
    def error_rate_threshold_used(self):
        '''

        self.aligned_segment.set_tag(tag = 'et', value = str(self.error_rate_threshold), value_type='Z', replace=True)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.error_rate_threshold_used
    
    @property
    def prior_probabilities_used(self):
        '''

        self.aligned_segment.set_tag(tag = 'et', value = str(self.error_rate_threshold), value_type='Z', replace=True)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.prior_probabilities_used
    
    @property
    def error_rate_average_het_sites_used(self):
        '''

        self.aligned_segment.set_tag(tag = 'er', value = str(PhaseSet_max.error_rate_average_het_sites_used), value_type='Z', replace=True)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.error_rate_average_het_sites
    
    @property
    def per_base_mismatch_rate_used(self):
        '''

        self.aligned_segment.set_tag(tag = 'er', value = str(PhaseSet_max.error_rate_average_het_sites_used), value_type='Z', replace=True)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.per_base_mismatch_rate
    
    @property
    def total_hets_analyzed_favored_haplotype(self):
        '''
        Returns the total number of heterozygous positions that were included in the phasing evaluation for this aligned_segment
        that resulted in the most favorable likelihood (ie: positions containing heterozygous SNVs that align to the
        aligned_segment of interest.)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.total_hets_analyzed
    
    @property
    def total_hets_favored_haplotype(self):
        '''
        Returns the total number of heterozygous positions that overlapped this aligned_segment before filtering the sites with
        '-'. This set of het sites was the most favorable phase set that resulted in the most favorable likelihood (
        ie: positions containing heterozygous SNVs that align to the aligned_segment of interest.)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.total_hets
    
    @property
    def total_matches_to_favored_haplotype(self):
        '''
        Returns the total number of matches to the favored haplotype sequence at heterozygous positions that
        overlapped this aligned_segment. This phased set of het sites resulted in the most favorable likelihood for this aligned_segment.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.matches[self.PhaseSet_max.max_phase - 1]
    
    @property
    def total_non_matches_to_favored_haplotype(self):
        '''
        Returns the total number of non-matches to the favored haplotype sequence at heterozygous positions that
        overlapped this aligned_segment. This phased set of het sites resulted in the most favorable likelihood for this aligned_segment.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.non_matches[self.PhaseSet_max.max_phase - 1]
    
    def prior_probability_for_haplotype_i(self, i):
        '''
        self.aligned_segment.set_tag(tag = 'pr', value = str(self.prior_probabilities), value_type='Z', replace=True)


        prior_probability_haplotype_i = P(haplotype_i)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.prior_probabilities[i - 1]
    
    @property
    def prior_probability_for_favored_haplotype(self):
        '''
        self.aligned_segment.set_tag(tag = 'pr', value = str(self.prior_probabilities), value_type='Z', replace=True)


        prior_probability_haplotype_i = P(haplotype_i)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.prior_probabilities[self.max_phase - 1]
    
    def calculate_log_likelihood_read_given_haplotype_i(self, i):
        """

        Args:
            i:

        Returns:

        """
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.log_probability_read_given_haplotype_i[i - 1]
    
    def calculate_likelihood_of_read_given_haplotype_i(self, i):
        """

        Args:
            i:

        Returns:

        """
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return 10**self.PhaseSet_max.log_probability_read_given_haplotype_i[i - 1]
    
    @property
    def calculate_log_likelihood_read_given_favored_haplotype(self):
        """

        Returns:

        """
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.PhaseSet_max.log_probability_read_given_haplotype_i[self.max_phase - 1]
    
    @property
    def calculate_likelihood_of_read_given_favored_haplotype(self):
        """

        Returns:

        """
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return 10**self.PhaseSet_max.log_probability_read_given_haplotype_i[self.max_phase - 1]
    
    @property
    def bayes_numerator_for_favored_haplotype(self):
        '''
        P(aligned_segment|haplotype_i)P(haplotype_i)
        Bayes numerator. likelihood * prior
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.calculate_likelihood_of_read_given_haplotype_i(
                    self.max_phase - 1
                    ) * self.prior_probability_for_haplotype_i(self.max_phase - 1)
        # if self.likelihood_for_haplotype_i(i) and self.prior_probability_for_haplotype_i(i):
    
    @property
    def log_likelihood_ratio_for_favored_haplotype_vs_not_favored_haplotype(self):
        '''

        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(
                    self.max_phase - 1
                    ):
                return math.log10(self.bayes_numerator_for_favored_haplotype) - math.log10(
                        sum(
                                [self.bayes_numerator_for_haplotype_i(_i - 1) for _i in
                                 range(1, self.ploidy_phase_set + 1) if
                                 _i - 1 != self.max_phase - 1]
                                )
                        )
    
    @property
    def posterior_probability_for_favored_haplotype(self):
        '''
        # posterior_probability_haplotype_i = P(haplotype_i|aligned_segment)

        # aligned_segment.set_tag(tag = 'fp', value = str(phased_read.posterior_probability_for_haplotype_i(c)), value_type='Z',
        replace=True) # aligned_segment.set_tag(tag = 'pp', value = str(phased_read.posterior_probability_for_haplotype_i(int(
        aligned_segment.get_tag('oa')))), value_type='Z', replace=True)


        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(
                    self.max_phase - 1
                    ):
                return self.bayes_numerator_for_haplotype_i(
                        self.max_phase - 1
                        ) / self.total_probability_of_read_given_haplotypes
    
    @property
    def log_posterior_probability_for_favored_haplotype(self):
        '''
        posterior_probability_haplotype_i = P(haplotype_i|aligned_segment)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(
                    self.max_phase - 1
                    ):
                return math.log10(self.bayes_numerator_for_haplotype_i(self.max_phase - 1)) - math.log10(
                        self.total_probability_of_read_given_haplotypes
                        )
    
    def bayes_numerator_for_haplotype_i(self, i):
        '''
        P(aligned_segment|haplotype_i)P(haplotype_i)
        Bayes numerator. likelihood * prior
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return self.calculate_likelihood_of_read_given_haplotype_i(i - 1) * self.prior_probability_for_haplotype_i(
                    i - 1
                    )
        # if self.likelihood_for_haplotype_i(i) and self.prior_probability_for_haplotype_i(i):
    
    def log_likelihood_ratio_for_haplotype_i_vs_not_haplotype_i(self, i):
        '''

        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(i - 1):
                return math.log10(self.bayes_numerator_for_haplotype_i(i - 1)) - math.log10(
                        sum(
                                [self.bayes_numerator_for_haplotype_i(_i - 1) for _i in
                                 range(1, self.ploidy_phase_set + 1) if
                                 _i - 1 != i - 1]
                                )
                        )
    
    def posterior_probability_for_haplotype_i(self, i):
        '''
        # posterior_probability_haplotype_i = P(haplotype_i|aligned_segment)

        # aligned_segment.set_tag(tag = 'fp', value = str(phased_read.posterior_probability_for_haplotype_i(c)), value_type='Z',
        replace=True) # aligned_segment.set_tag(tag = 'pp', value = str(phased_read.posterior_probability_for_haplotype_i(int(
        aligned_segment.get_tag('oa')))), value_type='Z', replace=True)


        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(i - 1):
                return self.bayes_numerator_for_haplotype_i(i - 1) / self.total_probability_of_read_given_haplotypes
    
    def log_posterior_probability_for_haplotype_i(self, i):
        '''
        posterior_probability_haplotype_i = P(haplotype_i|aligned_segment)
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            if self.total_probability_of_read_given_haplotypes and self.bayes_numerator_for_haplotype_i(i - 1):
                return math.log10(self.bayes_numerator_for_haplotype_i(i - 1)) - math.log10(
                        self.total_probability_of_read_given_haplotypes
                        )
    
    @property
    def total_probability_of_read_given_haplotypes(self):
        '''
        total_probability_of_read_given_haplotypes = sum(P(aligned_segment|haplotype_i)P(haplotype_i)) for i = 1,...,ploidy Our
        model assumes that the haplotypes listed in the VCF file are the ground truth molecular composition of the
        populations from which these reads were sampled from physically.
        '''
        if self.PhaseSet_max is None:
            return
        elif self.PhaseSet_max.phase == 'Nonphasable':
            return 'Nonphasable'
        else:
            return sum([self.bayes_numerator_for_haplotype_i(i) for i in range(1, self.ploidy_phase_set + 1)])

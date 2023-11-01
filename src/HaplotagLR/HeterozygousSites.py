from collections import defaultdict
from decimal import Decimal
from fractions import Fraction
from statistics import mean
from typing import Any, List, Optional, Tuple, Union

from Bio import Align


# noinspection PyUnusedClass
class HeterozygousSite:
    #probs_read: List[Union[float, Any]]
    #aligned_het_bases: defaultdict[Any, dict]
    #total_non_matches: int
    #aligned_het_bases: defaultdict[Any, dict]
    #local_realignment_counter: defaultdict[Any, int]
    # best_haplotype: List[Any]
    #alignments_favored_score_for_vary_lengths: List[Tuple[Optional[Any], List[Any]]]
    #average_score_diff_for_local_alignments: Union[Union[float, Decimal, Fraction], Any]
    sample: str
    vcf_record: object
    
    def __init__(self, vcf_record, sample: object, aligned_segment, liftover_converters = None):
        """
        Args:
            vcf_record:
            sample (object):
            aligned_segment:
            liftover_converters:
        """
        self.total_non_matches = None
        self.probs_read = None
        self.haplotype = None
        self.vcf_record = vcf_record
        self.sample = sample
        self.aligned_segment = aligned_segment
        if liftover_converters:
            self.liftover_converters = liftover_converters
        else:
            self.liftover_converters = None
    
    @property
    def is_snv(self) -> bool:
        """
        Returns:
            bool:
        """
        SNV: bool = True
        for haplotype in range(0, self.ploidy_het_site):
            if not len(self.alleles[haplotype]) == 1:
                SNV = False
        return SNV
    
    @property
    def is_match(self) -> bool:
        match = False
        for haplotype in range(0, self.ploidy_het_site):
            if self.alleles[haplotype] == self.read_base:
                match = True
        
        return match
    
    @property
    def is_match_to_phase(self) -> List[bool]:
        """
        Returns:
            object:
        """
        matches: List[bool] = []
        for haplotype in range(0, self.ploidy_het_site):
            match: bool = False
            if self.alleles[haplotype] == self.read_base:
                match = True
            matches.append(match)
        
        return matches
    
    @property
    def matches_single_phase(self) -> bool:
        """
        Returns:
            bool:
        """
        matches = []
        for haplotype in range(0, self.ploidy_het_site):
            match = False
            if self.alleles[haplotype] == self.read_base:
                match = True
            if match:
                matches.append(haplotype + 1)
        if len(matches) == 1:
            return True
        else:
            return False
    
    @property
    def haplotype_matches(self) -> object:
        """
        Returns:
            List[int]:
        """
        matches: List[int] = []
        haplotype: int
        for haplotype in range(0, self.ploidy_het_site):
            match = False
            if self.alleles[haplotype] == self.read_base:
                match = True
            if match:
                matches.append(haplotype + 1)
        return matches
    
    def is_non_match(self, haplotype: int = None) -> bool:
        """Args: int
            haplotype:

        Args:
            haplotype (int):

        Returns:
            bool:
        """
        
        # haplotype: int
        # for haplotype in range(0, self.ploidy_het_site):
        #     if not self.alleles[haplotype] == self.read_base:
        #         return False
        
        return self.alleles[haplotype] == self.read_base

    @property
    def ploidy_het_site(self) -> int:
    
        """INPUT ATTRIBUTE = FALSE OUTPUT ATTRIBUTE = TRUE :returns: :rtype:
        int
        """
    
        return self.vcf_record.samples[str(self.sample)].alleles
    
        # assert isinstance(self.alleles, object)
        # return len(self.alleles)
    
    @property
    def is_heterozygous(self) -> bool:
        """
        Returns:
            bool:
        """
        return min(self.vcf_record.samples[str(self.sample)].allele_indices) < max(
                self.vcf_record.samples[str(self.sample)].allele_indices
                )
    
    @property
    def reference_position(self) -> int:
        """
        Returns:
            int:
        """
        return self.vcf_record.pos - 1
    
    @property
    def reference_name(self) -> str:
        """
        Returns:
            str:
        """
        return str(self.aligned_segment.reference_name)
    
    @property
    def alleles(self) -> List[str]:
        """
        Returns:
            List[str]:
        """
        return self.vcf_record.samples[str(self.sample)].alleles
    
    @property
    def is_phased(self) -> bool:
        """
        Returns:
            bool:
        """
        return self.vcf_record.samples[str(self.sample)].phased
    
    @property
    def gapped_alignment_position(self) -> int:
        """
        Returns:
            int:
        """
        
        self._generate_gapped_alignment_position()
        
        return self._gapped_alignment_position
    
    @property
    def read_base(self) -> str:
        """
        Returns:
            str:
        """
        return self.aligned_segment.query_sequence[self.gapped_alignment_position]
    
    @property
    def read_base_quality(self) -> bytes:
        """
        Returns:
            bytes:
        """
        return self.aligned_segment.query_qualities[self.gapped_alignment_position]
    
    @property
    def read_base_quality_error_rate(self) -> float:
        """
        Returns:
            float:
        """
        return 10**(-0.1 * self.aligned_segment.query_qualities[self.gapped_alignment_position])
    
    # def read_base_quality_error_rate_for_offset(self, offset: object) -> object:
    #     """
    #
    #     Args:
    #         offset:
    #
    #     Returns:
    #
    #     """
    #     return 10**(-0.1 * self.aligned_segment.query_qualities[self.gapped_alignment_position + offset])
    
    # @property
    # def correct_phase(self) -> bool:
    #     """
    #
    #     SIMULATED SUBCLASS
    #     Returns:
    #         bool:
    #
    #     """
    #     if not self.aligned_segment.has_tag('oa'):
    #         return
    #     return self.aligned_segment.get_tag('oa')
    
    @property
    def max_allele_length(self) -> int:
        """
        Returns:
            int:
        """
        return max([len(allele) for allele in self.alleles])
    
    @property
    def allele_length(self) -> List[int]:
        return [len(allele) for allele in self.alleles]
    
    def allele_length_for_haplotype_i(self, haplotype) -> int:
        """
        Args:
            haplotype: int:

        Returns:
            int:
        """
        return [len(allele) for allele in self.alleles][int(haplotype) - 1]
    
    # @property
    def query_read_sequence(self, length: int) -> object:
        """
        Args:
            length (int):

        Returns:
            List[str]:
        """
        if self.gapped_alignment_position - length < 0:
            start = 0
        else:
            start = self.gapped_alignment_position - length
        
        if self.gapped_alignment_position + length - self.aligned_segment.query_length > 0:
            end = self.aligned_segment.query_length
        else:
            end = self.gapped_alignment_position + length
        
        return self.aligned_segment.query_sequence[start:end]
    
    # @property
    def target_ref_sequence(self, haplotype: object, length: object) -> object:
        """
        Args:
            haplotype (object): int
            length (object): int

        Returns:
            str:
        """
        if self.liftover_converters is None:
            return
        elif str(haplotype) in self.liftover_converters:
            # start = self.liftover_converters[str(haplotype)]['converter'].convert_coordinate(str(
            # self.reference_name), self.reference_position-length) print('start',start) if start: start = start[0][
            # 1] else: multiplier = 1 while len(start) == 0: start = self.liftover_converters[str(haplotype)][
            # 'converter'].convert_coordinate(str(self.reference_name), self.reference_position-(length*multiplier))
            # print('start',start, 'multiplier',multiplier) multiplier+=1 start = start[0][1] end =
            # self.liftover_converters[str(haplotype)]['converter'].convert_coordinate(str(self.reference_name),
            # self.reference_position+length) print('end',end) if end: end = end[0][1] else: multiplier = 1 while
            # len(end) == 0: end = self.liftover_converters[str(haplotype)]['converter'].convert_coordinate(str(
            # self.reference_name), self.reference_position+(length*multiplier)) print('end',end, 'multiplier',
            # multiplier) multiplier+=1 end = end[0][1] print('start',start) print('end',end)
            reference = str(self.reference_name)
            mid: object = self.liftover_converters[str(haplotype)]['converter'].convert_coordinate(
                    str(self.reference_name),
                    self.reference_position
                    )
            if len(mid) == 0:
                print(mid)
                multiplier = 1
                while len(mid) == 0:
                    mid = self.liftover_converters[str(haplotype)]['converter'].convert_coordinate(
                            str(self.reference_name), self.reference_position + multiplier
                            )
                    print('mid', mid, 'multiplier', multiplier)
                    multiplier += 1
            mid = mid[0][1]
            sequence: str = self.liftover_converters[str(haplotype)]['reference_fasta'].fetch(
                    reference = str(self.reference_name), start = mid - length, end = mid + length
                    )
            return sequence
    
    def local_realignment(self, length: object, upper_flag: object = True) -> object:
        """
        Args:
            length (object): int
            upper_flag (object): bool

        Returns:
            object:
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.target_open_gap_score = -1
        aligner.query_open_gap_score = -1
        aligner.mismatch_score = -0.5
        aligner.extend_gap_score = -0.5
        haplotype_favored = None
        query = self.query_read_sequence(length)
        
        alignments = defaultdict(list)
        for haplotype in range(1, self.ploidy_het_site + 1):
            target = self.target_ref_sequence(str(haplotype), length)
            if upper_flag:
                alignments[str(haplotype)] = aligner.align(str(target).upper(), str(query).upper())
            else:
                alignments[str(haplotype)] = aligner.align(target, query)
        
        return alignments
    
    # @property
    @property
    def calculate_average_score_diff_for_local_alignments(self) -> object:
        """
        Returns:
            float:
        """
        average_score_diff_for_local_alignments: Union[Union[float, Decimal, Fraction], Any]
        
        length: int
        self.average_score_diff_for_local_alignments = mean(
                [pair[1][0] for pair in [
                        calculate_favored_haplotype_from_local_alignment(self.local_realignment(length)) for length
                        in
                        range(5, 100, 5)] if len(pair[1]) > 0]
                )
        return self.average_score_diff_for_local_alignments
    
    # @property
    @property
    def automated_evaluation_of_local_alignments(self) -> object:
        """
        Returns:
            object:
        """
        alignments_favored_score_for_vary_lengths: List[Tuple[Optional[Any], List[Any]]]
        
        self.alignments_favored_score_for_vary_lengths = []
        for length in range(5, 100, 5):
            self.alignments_favored_score_for_vary_lengths.append(
                    calculate_favored_haplotype_from_local_alignment(self.local_realignment(length))
                    )
        self.sorted_pos = sorted(
                [pos for pos in zip(self.alignments_favored_score_for_vary_lengths, list(range(5, 100, 5))) if
                 pos[0][0]],
                key = lambda x:x[1], reverse = False
                )
        return self.alignments_favored_score_for_vary_lengths
    
    # @property
    @property
    def select_best_haplotype(self) -> int:
        """
        Returns:
            int:
        """
        best_haplotype: List[Any]
        self.automated_evaluation_of_local_alignments()
        self.count_best_haplotypes
        self.best_haplotype = []
        for haplotype in self.local_realignment_counter.items():
            haplotype: int
            if haplotype[1] == max(
                    [self.local_realignment_counter[haplotype] for haplotype in
                     self.local_realignment_counter]
                    ):
                self.best_haplotype.append(haplotype[0])
        if len(self.best_haplotype) == 1:
            return self.best_haplotype[0]
    
    # @property
    @property
    def count_best_haplotypes(self) -> int:
        """
        Returns:
            int:
        """
        local_realignment_counter: defaultdict[Any, int]
        self.local_realignment_counter = defaultdict(int)
        R: int
        for R in self.alignments_favored_score_for_vary_lengths:
            if 0 < len(R[1]):
                self.local_realignment_counter[R[0]] += R[1][0]
        if [self.local_realignment_counter[haplotype] for haplotype in self.local_realignment_counter if
            self.local_realignment_counter[haplotype]]:
            return max([self.local_realignment_counter[haplotype] for haplotype in self.local_realignment_counter])
    
    def calculate_aligned_het_bases_from_local_alignment(self, alignments: object) -> object:
        """
        Args:
            alignments (object): List[object]:

        Returns:
            int:
        """
        aligned_het_bases: defaultdict[Any, dict]
        self.aligned_het_bases: defaultdict[Any, dict] = defaultdict(dict)
        position: int
        for position in range(0, self.max_allele_length):
            query_position: int = self.gapped_alignment_position + position
            haplotype: object
            for haplotype in alignments:
                self.aligned_het_bases[str(haplotype)][query_position] = defaultdict(dict)
                
                # aligned_het_bases[haplotype][query_position]['target_base'] = defaultdict(list)
                # aligned_het_bases[haplotype][query_position]['query_base'] = defaultdict(list)
                
                index_length = int(len(alignments[str(haplotype)][0].target) / 2) + position
                if int(len(alignments[str(haplotype)][0].target) / 2) > position:
                    target_base = alignments[str(haplotype)][0].target[index_length]
                else:
                    target_base = '-'
                
                query_base: str = '-'
                for block in zip(alignments[str(haplotype)][0].aligned[0], alignments[str(haplotype)][0].aligned[1]):
                    if index_length in range(block[0][0], block[0][1]):
                        query_base = alignments[str(haplotype)][0].query[index_length - block[0][0] + block[1][0]]
                
                self.aligned_het_bases[haplotype][query_position][
                    'error_rate_quality_score'] = self.read_base_quality_error_rate_for_offset(position)
                self.aligned_het_bases[haplotype][query_position]['target_base'] = target_base
                self.aligned_het_bases[haplotype][query_position]['query_base'] = query_base
                self.aligned_het_bases[haplotype][query_position]['match'] = query_base == target_base
        
        return self.aligned_het_bases
    
    def calculate_probs(self, haplotype: object = None, force_recalc = False) -> float:
        """
        Args:
            haplotype (object):

        Returns:
            float:
        """
        # Use cached probs_read if possible.
        if len(self.probs_read) > 0 and not force_recalc:
            return sum(self.probs_read)
        
        total_non_matches: int = 0
        self.total_matches = 0
        #self.total_non_matches = 0
        aligned_het_bases: defaultdict[Any, dict]
        
        var = self.automated_evaluation_of_local_alignments
        if haplotype is None:
            haplotype = self.select_best_haplotype
        else:
            haplotype = haplotype
        
        var = self.automated_evaluation_of_local_alignments
        
        self.probs_read: List[Union[float, Any]] = []
        
        lengths: List[Any] = []
        for pos in self.sorted_pos:
            if pos[0][0] == str(haplotype):
                lengths.append(pos[1])
        if len(lengths) > 0:
            aligned_het_bases = self.calculate_aligned_het_bases_from_local_alignment(
                    self.local_realignment(lengths[0])
                    )
            self.aligned_het_bases = aligned_het_bases
            for offset in range(0, int(self.allele_length_for_haplotype_i(haplotype))):
                if aligned_het_bases[str(haplotype)][self.gapped_alignment_position + offset]['match']:
                    self.probs_read.append(
                            aligned_het_bases[str(haplotype)][self.gapped_alignment_position + offset][
                                'error_rate_quality_score']
                            )
                    self.total_matches += 1
                elif not aligned_het_bases[str(haplotype)][self.gapped_alignment_position + offset]['match']:
                    self.probs_read.append(
                            float(
                                    aligned_het_bases[str(haplotype)][self.gapped_alignment_position + offset][
                                        'error_rate_quality_score']
                                    ) / 3
                            )
                    self.total_non_matches += 1
            
            return sum(self.probs_read)
    
    @property
    def _generate_gapped_alignment_position(self) -> object:
    
        # reference_start, cigartuples, reference_positions):
        """
        # self.gapped_alignment_positions =
        self._generate_gapped_alignment(self.aligned_segment.reference_start,
        self.aligned_segment.cigartuples, self.reference_positions)

        Returns:
            object:
        """
    
        gapped_alignment_positions = []
        reference_start = self.aligned_segment.reference_start
        cigartuples = self.aligned_segment.cigartuples
        reference_positions: List[int] = [self.reference_position]
        
        i = reference_start  # self.aligned_segment.reference_start
        j = 0
        h = 0
        if not cigartuples:
            return gapped_alignment_positions
        
        cigartuple: tuple
        for cigartuple in cigartuples:  # self.aligned_segment.cigartuples:
            j1 = j
            if h >= len(reference_positions):
                break
            elif cigartuple[0] in [0, 7, 8]:
                i += cigartuple[1]
                j += cigartuple[1]
            elif cigartuple[0] in [1, 4]:
                j += cigartuple[1]
            elif cigartuple[0] in [2, 3]:
                i += cigartuple[1]
            while i > reference_positions[h]:
                if cigartuple[0] in [0, 7, 8]:
                    gapped_alignment_positions.append(j1 + (cigartuple[1] - (i - reference_positions[h])))
                    if h == len(reference_positions) - 1:
                        h += 1
                        break
                elif cigartuple[0] == 2:
                    gapped_alignment_positions.append(j1 + (cigartuple[1] - (i - reference_positions[h])))
                    if h == len(reference_positions) - 1:
                        h += 1
                        break
                if h >= len(reference_positions) - 1:
                    break
                else:
                    h += 1

        if len(gapped_alignment_positions) > 0:
            self._gapped_alignment_position = gapped_alignment_positions[0]
        return gapped_alignment_positions
    
    @property
    def correct_phase(self) -> object:
        """SIMULATED SUBCLASS :returns: :rtype: int"""
        # if not self.aligned_segment.has_tag('oa'):
        #     return
        return self.haplotype
    
    def read_base_quality_error_rate_for_offset(self, position):
        """
        Args:
            position:
        """
        pass


def calculate_favored_haplotype_from_local_alignment(alignments: List[object]) -> int:
    """
    Args:
        alignments (object):
    """
    
    diff: List[Any] = []
    scores: defaultdict[Any, list] = defaultdict(list)
    haplotype: object
    score: object
    for haplotype in alignments:
        scores[alignments[haplotype].score].append(haplotype)
    if len(scores[max(scores.keys())]) == 1:
        haplotype_favored: int = scores[max(scores.keys())][0]
        diff = []
        
        for score in scores.keys():
            if score != max(scores.keys()):
                diff.append(max(scores.keys()) - score)
    
    return haplotype_favored

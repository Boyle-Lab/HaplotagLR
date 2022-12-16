# coding=utf-8
import sys
import math
from collections import defaultdict
from statistics import mean
from typing import Any, List
import pysam
import numpy as np

import LRphase.HeterozygousSites as HeterozygousSite

# Minimum floating-point value > 0
# See http://www.seunosewa.com/2009/06/floating-point-min-max.html
minfloat = 4.9406564584124654e-324
log10minfloat = math.log10(minfloat)

def powlaw_modified(x: object, a: object = 4.5, xmin: object = 2) -> object:
    """

    Args:
        x:
        a:
        xmin:

    Returns:

    """
    return ((a - 1) / xmin) * math.pow(x / xmin, -1 * a)


"""
def _multcoeff(*args):
    '''
    # self.result = _multcoeff(*args)
    Return the multinomial coefficient
    (n1 + n2 + ...)! / n1! / n2! / ...
    '''
    if not args:  # no parameters
        return 1
    # Find and store the index of the largest parameter so we can skip
    #   it (for efficiency)
    skipndx = args.index(max(args))
    newargs = args[:skipndx] + args[skipndx + 1:]
    
    result = 1
    num = args[skipndx] + 1  # a factor in the numerator
    for n in newargs:
        for den in range(1, n + 1):  # a factor in the denominator
            result = result * num // den
            num += 1
    return result
"""

# Cache multinomial coefficients for efficiency
_multcoeff_cache = {}

def _multcoeff(n_match, n_mismatch):
    """
    Calculates multinomial coefficient according to the formula
    n! / ((n-m)! * m!)  (As described in the LRphase article)
    where n = the total number of heterozygous sites in a region 
    (n_match + n_mismatch) and m = the number of heterozygous sites
    matching the given haplotype (n_match). Values are cached to
    improve performance.
    """
    if (n_match, n_mismatch) not in _multcoeff_cache.keys():
        _multcoeff_cache[(n_match, n_mismatch)] = math.log10(np.math.factorial(n_match+n_mismatch) / (np.math.factorial(n_mismatch) * np.math.factorial(n_match)))
    return _multcoeff_cache[(n_match, n_mismatch)]
    

def _read_base_error_rates_from_quality_scores(gapped_alignment_positions, query_qualities):
    """# self.read_base_error_rates, self.error_rate_average_het_sites =
    self._read_base_error_rates_from_quality_scores(self.gapped_alignment_positions,
    self.aligned_segment.query_qualities) """
    
    ############ Extract base qualities for heterozygous sites in aligned_segment  ############
    
    read_base_error_rates = []
    error_rate_total_het_sites = 0
    total_het_sites_count = 0
    
    for position in gapped_alignment_positions:
        if isinstance(position, str):
            if position[0] == '-':
                read_base_error_rates.append(1)
        else:
            er = 10**(-0.1 * query_qualities[position])
            read_base_error_rates.append(er)
            error_rate_total_het_sites += er
            total_het_sites_count += 1
    
    if total_het_sites_count > 0:
        error_rate_average_het_sites = error_rate_total_het_sites / total_het_sites_count
    else:
        error_rate_average_het_sites = 0
    return read_base_error_rates, error_rate_average_het_sites


def _read_bases_from_gapped_alignment(gapped_alignment_positions, query_sequence):
    """# self.read_bases = self._read_bases_from_gapped_alignment(self.gapped_alignment_positions,
    self.aligned_segment.query_sequence) """
    ############ Extract bases for heterozygous sites in aligned_segment  ############
    read_bases = []
    if query_sequence is None:
        return read_bases
    for position in gapped_alignment_positions:
        if isinstance(position, str):
            if position[0] == '-':
                read_bases.append('-')
        else:
            read_bases.append(query_sequence[position])
    return read_bases


def _read_base_error_rates_from_gapped_alignment(gapped_alignment_positions: object, per_base_error_rate):
    """# self.read_base_error_rates = self._read_base_error_rates_from_gapped_alignment(
    self.gapped_alignment_positions, self.per_base_error_rate) """
    read_base_error_rates = []
    # PhaseSet._estimate_per_base_error_rate()
    for position in gapped_alignment_positions:
        if isinstance(position, str):
            if position[0] == '-':
                read_base_error_rates.append(1)  # Should this be 0.5? Right now does not matter because we don't use gap positions in scoring.
        else:
            read_base_error_rates.append(per_base_error_rate)
    return read_base_error_rates


def _calculate_per_base_mismatch_rate(aligned_segment):
    """# self.per_base_mismatch_rate = self._calculate_per_base_mismatch_rate(self.aligned_segment)"""
    matches = 0
    i = 0
    for position in aligned_segment.get_aligned_pairs(with_seq = True):
        if position[1] is None:
            continue
        elif position[0] is None:
            continue
        else:
            if str(aligned_segment.query_sequence[position[0]]) == str(position[2]).upper():
                matches += 1
        i += 1
    per_base_mismatch_rate = float(1 - (matches / i))
    return per_base_mismatch_rate


def _calculate_local_per_base_mismatch_rate(aligned_segment, gapped_alignment_positions) -> object:
    """# self.read_base_error_rates, self.per_base_mismatch_rate, self.error_rate_average_het_sites =
    self._calculate_local_per_base_mismatch_rate(self.aligned_segment, self.gapped_alignment_positions)

    Args:
        aligned_segment:
        gapped_alignment_positions: """
    
    # read_base_error_rates = defaultdict(int)
    error_rate_total_het_sites = 0
    total_het_sites_count = 0
    matches = 0
    matches_alignment_query: List[int] = []
    read_base_error_rates = []
    # gapped_alignment_positions_error_rate = []
    i = 0
    for position in aligned_segment.get_aligned_pairs(with_seq = True):
        if position[1] is None:
            matches_alignment_query.append(0)
            continue
        elif position[0] is None:
            continue
        else:
            if str(aligned_segment.query_sequence[position[0]]) == str(position[2]).upper():
                matches += 1
                matches_alignment_query.append(1)
            else:
                matches_alignment_query.append(0)
        i += 1
    per_base_mismatch_rate = float(1 - (matches / i))
    
    for site in gapped_alignment_positions:
        if isinstance(site, str):
            if site[0] == '-':
                read_base_error_rates.append(1)
        else:
            try:
                read_base_error_rates.append(
                        round(
                                (2 * round(
                                        1 - (
                                                sum(matches_alignment_query[int(site) - 10:int(site) + 10]) / len(
                                                matches_alignment_query[int(site) - 10:int(site) + 10]
                                                )), 4
                                        ) + round(
                                        1 - (
                                                sum(matches_alignment_query[int(site) - 25:int(site) + 25]) / len(
                                                matches_alignment_query[int(site) - 25:int(site) + 25]
                                                )), 4
                                        ) + round(
                                        1 - (
                                                sum(matches_alignment_query[int(site) - 50:int(site) + 50]) / len(
                                                matches_alignment_query[int(site) - 50:int(site) + 50]
                                                )), 4
                                        )) / 4, 4
                                )
                        )
                error_rate_total_het_sites += round(
                        (2 * round(
                                1 - (
                                        sum(matches_alignment_query[int(site) - 10:int(site) + 10]) / len(
                                        matches_alignment_query[int(site) - 10:int(site) + 10]
                                        )), 4
                                ) + round(
                                1 - (
                                        sum(matches_alignment_query[int(site) - 25:int(site) + 25]) / len(
                                        matches_alignment_query[int(site) - 25:int(site) + 25]
                                        )), 4
                                ) + round(
                                1 - (
                                        sum(matches_alignment_query[int(site) - 50:int(site) + 50]) / len(
                                        matches_alignment_query[int(site) - 50:int(site) + 50]
                                        )), 4
                                )) / 4, 4
                        )
            except ZeroDivisionError:
                read_base_error_rates.append(per_base_mismatch_rate)
                error_rate_total_het_sites += per_base_mismatch_rate
            total_het_sites_count += 1
    
    if total_het_sites_count > 0:
        error_rate_average_het_sites = error_rate_total_het_sites / total_het_sites_count
    else:
        error_rate_average_het_sites = 0
    
    return read_base_error_rates, per_base_mismatch_rate, error_rate_average_het_sites


def _count_matches(
        read_bases: object, phased_alleles: object,
        read_base_error_rates: object = None, per_base_error_rate: float = None,
        multinomial_correction: object = None, pseudocount: float = 0.001
        ) -> object:
    """
    # self.log_probability_read_given_haplotype_i, self.non_matches, self.matches, self.total_hets,
    self.total_hets_analyzed, self.ploidy_phase_set = self._count_matches(self.read_bases, self.phased_alleles,
    self.read_base_error_rates, self.error_model)

    Args:
        read_bases:
        phased_alleles:
        read_base_error_rates:
        multinomial_correction:

    Returns:
        object: 
    """
    log_probability_read_given_haplotype_i: defaultdict[Any, int] = defaultdict(int)  # log10(P(data|phase_i))
    non_matches = defaultdict(int)
    matches = defaultdict(int)
    total_hets = 0
    total_hets_analyzed = 0
    ploidy_list = []
    ploidy_phase_set = 0

    # Sanity check for aberant invocation
    if per_base_error_rate is None and read_base_error_rates is None:
        # Throw a warning and return an empty result.
        sys.stderr.write("WARNING: _count_matches invoked with no values for per_base_error_rate or read_base_error_rates. Affected read(s) will be labelled nonphasable.")
        return log_probability_read_given_haplotype_i, non_matches, matches, total_hets, total_hets_analyzed, ploidy_phase_set

    if per_base_error_rate is not None:
        q = per_base_error_rate
        # We've seen weirdness with log(0) errors in this portion of the code, so we
        # need to explicitly check for q == 0.
        #if q == 0:
        # It turns out the errors are cropping up due to negative values being supplied
        # for per_base_error_rate. Don't know where these are coming from right now, or why,
        # but for now we can safeguard against this by checking for values < 0 too.
        if q <= 0:  # Debugging 
            q = sys.float_info.min
        match_log_prob = math.log10(1 - q)
        # Even with the above check, we've still seen log(0) errors here (which makes
        # no sense). Until we can figure out why we're still seeing q == 0, wrapping
        # this in a try/except block should prevent crashes.
        try:
            error_log_prob = math.log10(q / 3)
        except:
            sys.stderr.write("_count_matches: Math domain error encountered. Args: %s, %s, %s, %s, %s, %s\n" % (read_bases, phased_alleles, read_base_error_rates, per_base_error_rate, multinomial_correction, pseudocount))
            # For debugging purposes, return and flag one of the values with something we never expect to see otherwise:
            return 9999999999999999999999999, non_matches, matches, total_hets, total_hets_analyzed, ploidy_phase_set
            error_log_prob = math.log10(sys.float_info.min)

    # print(read_bases, phased_alleles, read_base_error_rates, error_model)
    for i in range(len(read_bases)):
        read_base = read_bases[i]
        ploidy = len(phased_alleles[i])
        ploidy_list.append(ploidy)
        
        if not read_base == '-': # Skip gap positions
            """
            The loop below does not work consistently. In its current state, it is inconsistent
            in how complex variants are handled (i.e., those where one or more alleles is >1bp)
            in that these are counted as mismatches if the allele for the present haplotype is
            a single bp, but are skipped entirely when the current allele has len > 1. This means
            that log-likelihoods are calculated using different numbers of heterozygous sites per
            haplotype, which is incorrect. This commit adds functionality to filter complex variants
            out of the variant list prior to even running this step. This is done by default when the
            PhaseSet object is initialized. Given this change in behavior, the current implementation
            of this loop produces correct and consistent scores. However, this should be modified in
            future releases to calculate correct scores given all types of potential variants. In
            short, I'm leaving this as-is for now and relying on only simple SNPs being passed in
            in phased_alleles.
            """
            if read_base_error_rates is not None:
                """
                This scoring mode equates to a generalized multinomial model, where
                the error probability is allowed to vary from site to site. Therefore,
                the current implementation of _multcoeff will not produce a valid
                coefficient (see https://fractional-calculus.com/multinomial_theorem.pdf)
                Before we can use these values, we need to generalize the multinomial
                coefficient as described in the referenced pdf. Briefly, _multcoeff needs
                to calculate fractional coefficients for each observed rate and sum the
                values.                
                """
                # Precomputing these does not save time because it still requires
                # a loop and these are only used once.
                q = read_base_error_rates[i] + pseudocount
                if q <= 0:
                    q = sys.float_info.min
                match_log_prob = math.log10(1-q)
                try:
                    error_log_prob = math.log10(q/3)
                except:
                    error_log_prob = math.log10(sys.float_info.min)
                ## TO-DO: calculate read-specific generalized multinomial coefficient
                
            for haplotype in range(0, ploidy):
                if len(phased_alleles[i][haplotype]) > 1:  # Probably to exclude indels/rearrangements
                    continue
                total_hets_analyzed += 0.5  # These are counted twice. Should prob be 1/ploidy instead of hard-coded 0.5.
                if read_base == phased_alleles[i][haplotype]:
                    log_probability_read_given_haplotype_i[haplotype] += match_log_prob
                    matches[haplotype] += 1
                else:
                    log_probability_read_given_haplotype_i[haplotype] += error_log_prob
                    non_matches[haplotype] += 1
        total_hets += 1
    """
    sys.stderr.write("Total hets: %d, total_analyzed: %d, matches[0]: %d, mismatches[0]: %d, matches[1], %d, mismatches[1]: %d, q: %.5f\n" %
                     (total_hets, total_hets_analyzed,
                      matches[0], non_matches[0],
                      matches[1], non_matches[1], q)
    )
    """
    if total_hets_analyzed > 0:
        ploidy_phase_set = max(ploidy_list)
    
    if multinomial_correction:
        haplotype: int
        for haplotype in range(0, ploidy_phase_set):
            multinomial = _multcoeff(int(matches[haplotype]), int(non_matches[haplotype]))
            log_probability_read_given_haplotype_i[haplotype] = log_probability_read_given_haplotype_i[
                                                                    haplotype] + multinomial
            #sys.stderr.write("Haplotype %d log_prob: %.5f, prior: %.5f\n" % (haplotype, log_probability_read_given_haplotype_i[haplotype], multinomial))
    return log_probability_read_given_haplotype_i, non_matches, matches, total_hets, total_hets_analyzed, ploidy_phase_set


def _estimate_per_base_error_rate(aligned_segment: pysam.AlignedSegment) -> object:
    """# self.per_base_error_rate, self.error_rate_calculation = self._estimate_per_base_error_rate(
    self.aligned_segment) """
    error_rate_calculation = []
    try:
        per_base_error_rate = aligned_segment.get_tag(
                tag = 'de'
                )  # Gap-compressed per-base sequence divergence from minimap2 output
        error_rate_calculation.append('de tag')
    except LookupError:
        error_rate_calculation.append('no de tag')
        try:
            per_base_error_rate = aligned_segment.get_tag(tag = 'NM') / aligned_segment.query_alignment_length
            error_rate_calculation.append('NM tag')
        except LookupError:
            error_rate_calculation.append('no NM tag')
            try:
                per_base_error_rate = _calculate_per_base_mismatch_rate(aligned_segment)
                error_rate_calculation.append('_calculate_per_base_mismatch_rate success')
            except LookupError:
                error_rate_calculation.append('_calculate_per_base_mismatch_rate failed')
    if per_base_error_rate == 0:
        per_base_error_rate = 0.001
    return per_base_error_rate, error_rate_calculation


def _get_phased_alleles(vcf_records, sample, use_complex_variants = True, use_gap_variants = True):
    '''self.phased_alleles = self._get_phased_alleles(self.vcf_records, self.sample) '''
    provisional_recs = [vcf_record.samples[sample].alleles for vcf_record in vcf_records]

    if use_complex_variants and use_gap_variants:
        return provisional_recs
    elif len(provisional_recs) == 0:
        return provisional_recs

    haplotypes = range(len(provisional_recs[0]))
    ret = []
    for rec in provisional_recs:
        use_rec = True
        for haplotype in haplotypes:
            if not use_complex_variants and len(rec[haplotype]) > 1: # Drop complex variants
                use_rec = False
                break
            if not use_gap_variants and rec[haplotype] == '-': # Drop gap variants
                use_rec = False
                break
        if use_rec:
            ret.append(rec)

    return ret


def _get_reference_positions(vcf_records, sample,  use_complex_variants = True, use_gap_variants = True):
    '''self.reference_positions = self._get_reference_positions(self.vcf_records) '''
    provisional_pos = [vcf_record.pos - 1 for vcf_record in vcf_records] # Converts pos to zero-based coord. frame
    provisional_recs = [vcf_record.samples[sample].alleles for vcf_record in vcf_records]

    if use_complex_variants and use_gap_variants:
        return provisional_pos
    elif len(provisional_recs) == 0:
        return provisional_pos

    haplotypes = range(len(provisional_recs[0]))
    ret = []
    for i in range(len(provisional_recs)):
        rec = provisional_recs[i]
        use_rec = True
        for haplotype in haplotypes:
            if not use_complex_variants and len(rec[haplotype]) > 1:
                use_rec = False
                break
            if not use_gap_variants and rec[haplotype] == '-':
                use_rec = False
                break
        if use_rec:
            ret.append(provisional_pos[i])

    return ret


def _generate_gapped_alignment(reference_start, cigartuples, reference_positions) -> object:
    """# self.gapped_alignment_positions = self._generate_gapped_alignment(self.aligned_segment.reference_start,
    self.aligned_segment.cigartuples, self.reference_positions)

    Args:
        reference_start:
        cigartuples:
        reference_positions: """
    i: object = reference_start  # self.aligned_segment.reference_start
    j = 0
    h = 0
    gapped_alignment_positions = []
    
    if not cigartuples:
        return gapped_alignment_positions
    
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
                gapped_alignment_positions.append('-' + str(j1 + (cigartuple[1] - (i - reference_positions[h]))))
                if h == len(reference_positions) - 1:
                    h += 1
                    break
            if h >= len(reference_positions) - 1:
                break
            else:
                h += 1
    return gapped_alignment_positions


def _estimate_prior_probabilities(ploidy_phase_set = 2, prior_probabilities: object = None):
    """
    Returns log10 transformed priors.
    """
    if prior_probabilities is None:
        if ploidy_phase_set > 0:
            prior_probabilities = list(map(math.log10, [1 / ploidy_phase_set] * ploidy_phase_set))  # P(phase_i)
    return prior_probabilities


class PhaseSet:
    """

    """
    
    def __init__(
            self, phase_set, aligned_segment, vcf_records, sample, error_model = 0, error_rate_threshold = 0.01,
            multinomial_correction = True, prior_probabilities = None, liftover_converters = None, phase = 'Nonphasable',
            powlaw_alpha = 4.5, powlaw_xmin = 2, use_complex_variants = False, use_gap_variants = False
            ):
        self.aligned_segment = aligned_segment
        self.phase_set = phase_set
        self.vcf_records = vcf_records
        if liftover_converters:
            self.liftover_converters = liftover_converters
        else:
            self.liftover_converters = None
        # self.vcf_records
        self.sample = sample
        self.error_model = error_model
        self.error_rate_threshold = error_rate_threshold
        self.prior_probabilities = _estimate_prior_probabilities()  # Log10 transformed priors
        self.multinomial_correction = multinomial_correction
        self.gapped_alignment_positions = []
        self.read_bases = []
        self.phased_alleles = []
        self.reference_positions = []
        self.max_log_likelihood_ratio = None
        self.max_phase = None
        self.phase = phase
        self.powlaw_alpha = powlaw_alpha
        self.powlaw_xmin = powlaw_xmin
        
        self.use_complex_variants = use_complex_variants
        self.use_gap_variants = use_gap_variants
        
        self.total_hets = 0
        if not self.aligned_segment.is_unmapped: # or self.aligned_segment.is_aligned_to_contig_not_in_vcf:
            # Don't think second test is necessary because it should never yield len(self.vcf_records) > 0.
            # Current implementation of PhasedRead object sometimes leaves this attribute unset, causing
            # exceptions, so better if we can avoid it.
            if len(self.vcf_records) > 0:
                self.reference_positions = _get_reference_positions(self.vcf_records, self.sample, self.use_complex_variants, self.use_gap_variants)
                self.phased_alleles = _get_phased_alleles(self.vcf_records, self.sample, self.use_complex_variants, self.use_gap_variants)
                self.gapped_alignment_positions = _generate_gapped_alignment(
                        self.aligned_segment.reference_start,
                        self.aligned_segment.cigartuples,
                        self.reference_positions
                        )
                self.read_bases = _read_bases_from_gapped_alignment(
                        self.gapped_alignment_positions,
                        self.aligned_segment.query_sequence
                        )
        self.read_base_error_rates = None
        self.per_base_mismatch_rate = None
        self.error_rate_average_het_sites = None

        
    def __iter__(self):
        self.HeterozygousSites = []
        self.HeterozygousSites_processed_count = 0
        # self.alignment_files_read_counts = [bam_file.mapped+bam_file.unmapped for bam_file in self.bam_files]
        return self
    
    def __next__(self):
        if self.HeterozygousSites_processed_count < len(self.vcf_records):
            Heterozygous_Site = HeterozygousSite(
                    self.vcf_records[self.HeterozygousSites_processed_count], self.sample,
                    self.aligned_segment, self.liftover_converters
                    )
            self.HeterozygousSites.append(Heterozygous_Site)
            self.HeterozygousSites_processed_count += 1
            return Heterozygous_Site
        else:
            raise StopIteration()
    
    def solve_phase(
            self, error_model: object = 0, error_rate_threshold: object = 0.01, prior_probabilities: object = None,
            multinomial_correction: object = True
            ) -> object:
        """

        Args:
            error_model:
            error_rate_threshold:
            prior_probabilities:
            multinomial_correction:

        Returns:

        """
        if prior_probabilities is None:
            prior_probabilities = self.prior_probabilities
        if not self.aligned_segment.is_unmapped or self.aligned_segment.is_aligned_to_contig_not_in_vcf:
            if len(self.gapped_alignment_positions) > 0:

                """
                Note that, as of 12/15/2022, error_models > 0 are disabled in cli.py since they will not
                yield valid multinomial probabilities under the current implementation of _multcoeff.
                See further notes in _count_matches.
                """
                
                if error_model == 0:
                    self.per_base_mismatch_rate = self.aligned_segment.get_tag(tag = 'de')
                    self.error_rate_average_het_sites = self.aligned_segment.get_tag(tag = 'de')

                elif error_model < 5:
                    self.read_base_error_rates, self.error_rate_average_het_sites, self.per_base_mismatch_rate = self._retrieve_read_base_error_rates(
                        self.aligned_segment,
                        self.gapped_alignment_positions,
                        error_model
                    )

                if error_model < 5:
                    self.log_probability_read_given_haplotype_i, self.non_matches, self.matches, self.total_hets, self.total_hets_analyzed, self.ploidy_phase_set = _count_matches(
                        self.read_bases,
                        self.phased_alleles,
                        self.read_base_error_rates,
                        self.per_base_mismatch_rate,
                        multinomial_correction
                    )
                    # This is to watch for debug info from a log10(0) error:
                    if self.log_probability_read_given_haplotype_i == 9999999999999999999999999:
                        sys.stderr.write("%s\n" % (self.aligned_segment))
                    
                    self.log_likelihood_ratios, self.phase, self.max_phase, self.max_log_likelihood_ratio = _phasing_decision(
                        self.log_probability_read_given_haplotype_i,
                        self.ploidy_phase_set,
                        prior_probabilities,
                        error_rate_threshold,
                        self.total_hets_analyzed,
                        powlaw_alpha=self.powlaw_alpha,
                        powlaw_xmin=self.powlaw_xmin
                    )

                    #sys.stderr.write("het sites: %s, matchesm: %s, non_matches: %s\n" % (self.total_hets, self.matches, self.non_matches))
                
                elif error_model == 5:
                    haplotypes = {}
                    self.matches = 0
                    self.non_matches = 0
                    hets = [het for het in self]
                    self.total_hets = len(hets)
                    self.log_likelihood_ratios = []
                    log_likelihood_ratio = {}
                    self.ploidy_phase_set = int(mean([int(het.ploidy_het_site) for het in self]))
                    for haplotype in range(1, self.ploidy_phase_set + 1):
                        haplotypes[haplotype] = sum(
                                [het.calculate_probs(haplotype = haplotype) for het in hets if
                                 het.calculate_probs(haplotype = haplotype)] # Second call to calculate_probs uses cached values
                                )
                        # if len(haplotypes[haplotype])>0:
                        # print(haplotypes[haplotype])
                    # print('1', haplotypes)
                    # print('2', haplotypes.values())
                    # print('3',list(haplotypes.values()))
                    
                    for count, haplotype in enumerate(haplotypes):
                        self.log_likelihood_ratios.append(haplotypes[haplotype] - haplotypes[len(haplotypes) - count])
                    self.max_log_likelihood_ratio = max(self.log_likelihood_ratios)
                    self.max_phase = self.log_likelihood_ratios.index(self.max_log_likelihood_ratio) + 1
                    if self.max_log_likelihood_ratio > self.error_rate_threshold:
                        self.phase = self.max_phase
                    else:
                        self.phase = 'Unphased'
                    self.total_hets = len(hets)
                    self.total_hets_analyzed = len(hets)
                    
                    self.log_probability_read_given_haplotype_i = haplotypes
                    for het in hets:
                        het.calculate_probs()
                        if het.total_matches > 0:
                            self.matches += 1
                        elif het.total_non_matches > 0:
                            self.non_matches += 1
                    
                    # self.matches = sum([1 for het in self if het.calculate_probs()[str(self.max_phase)][
                    # het.gapped_alignment_position]['match']]) self.non_matches = sum([1 for het in self if not
                    # het.calculate_probs()[str(self.max_phase)][het.gapped_alignment_position]['match']])
                    # per_base_mismatch_rate = aligned_segment.get_tag(tag='NM')/aligned_segment.query_alignment_length
                    # read_base_error_rates = [het.read_base_quality_error_rate for het in
                    # self]self._read_base_error_rates_from_gapped_alignment(gapped_alignment_positions,
                    # per_base_mismatch_rate)
        
        self.error_model_used = error_model
        self.error_rate_threshold_used = error_rate_threshold
        self.prior_probabilities_used = prior_probabilities
        
        return  # self.max_log_likelihood_ratio, self.phase_set
    
    def _retrieve_read_base_error_rates(
            self, aligned_segment: object, gapped_alignment_positions: object,
            error_model: object
            ) -> object:
        """

        Args:
            aligned_segment:
            gapped_alignment_positions:
            error_model:

        Returns:

        """
        error_rate_average_het_sites = None
        read_base_error_rates = []
        per_base_mismatch_rate = None
        error_rate_average_het_sites = None
        
        if not self.aligned_segment.is_unmapped or self.aligned_segment.is_aligned_to_contig_not_in_vcf:
            if len(self.gapped_alignment_positions) > 0:

                if error_model == 0:
                    per_base_mismatch_rate = aligned_segment.get_tag(tag = 'de')
                    read_base_error_rates = _read_base_error_rates_from_gapped_alignment(
                        gapped_alignment_positions, per_base_mismatch_rate
                    )
                
                if error_model == 1:
                    read_base_error_rates, per_base_mismatch_rate, error_rate_average_het_sites = _calculate_local_per_base_mismatch_rate(
                        aligned_segment, gapped_alignment_positions
                    )
                    
                elif error_model == 2:
                    read_base_error_rates, error_rate_average_het_sites = _read_base_error_rates_from_quality_scores(
                        gapped_alignment_positions, aligned_segment.query_qualities
                    )
                    per_base_mismatch_rate = aligned_segment.get_tag(
                        tag = 'NM'
                    ) / aligned_segment.query_alignment_length
                    
                elif error_model == 3:
                    per_base_mismatch_rate = aligned_segment.get_tag(
                        tag = 'NM'
                    ) / aligned_segment.query_alignment_length
                    read_base_error_rates = _read_base_error_rates_from_gapped_alignment(
                        gapped_alignment_positions, per_base_mismatch_rate
                    )
                    
                elif error_model == 4:
                    per_base_mismatch_rate = _calculate_per_base_mismatch_rate(aligned_segment)
                    read_base_error_rates = _read_base_error_rates_from_gapped_alignment(
                        gapped_alignment_positions, per_base_mismatch_rate
                    )
                    
                if len(read_base_error_rates) > 0:
                    if len([error_rate for error_rate in read_base_error_rates if error_rate != 1]) > 0:
                        error_rate_average_het_sites = mean(
                                [error_rate for error_rate in read_base_error_rates if error_rate != 1]
                                )
        
        return read_base_error_rates, error_rate_average_het_sites, per_base_mismatch_rate


def _phasing_decision(
        log_probability_read_given_haplotype_i, ploidy_phase_set, prior_probabilities,
        error_rate_threshold, total_hets_analyzed, powlaw_alpha = 4.5, powlaw_xmin = 2,
        minfloat = log10minfloat
        ) -> object:
    """# self.log_likelihood_ratios, self.phase, self.max_phase, self.max_log_likelihood_ratio =
    self._phasing_decision(self.log_probability_read_given_haplotype_i, self.ploidy_phase_set,
    self.prior_probabilities, self.error_rate_threshold)

    Args:
        log_probability_read_given_haplotype_i:
        ploidy_phase_set:
        prior_probabilities:
        error_rate_threshold:
        total_hets_analyzed: """
    
    if total_hets_analyzed == 0:
        phase = 'Nonphasable'
        log_likelihood_ratios = None
        max_log_likelihood_ratio: float = None
        max_phase = 'Nonphasable'
        return log_likelihood_ratios, phase, max_phase, max_log_likelihood_ratio
    
    else:
        log_likelihood_ratios = []
        for haplotype in log_probability_read_given_haplotype_i:
            denominator = 0
            numerator = 0
            for j in range(0, ploidy_phase_set):
                if haplotype == j:
                    # Assumes priors are log-transformed! 
                    numerator += (log_probability_read_given_haplotype_i[haplotype] +
                                  prior_probabilities[haplotype])
                    """
                    try:
                        # Assumes priors are log-transformed!
                        numerator += (log_probability_read_given_haplotype_i[haplotype] +
                                      prior_probabilities[haplotype])
                    except ValueError:
                        numerator += log_probability_read_given_haplotype_i[haplotype]
                    """
                else:
                    # Assumes priors are log-transformed! 
                    denominator += 10**(log_probability_read_given_haplotype_i[j] + prior_probabilities[j])
            #sys.stderr.write("hap: %d, numer: %f, denom %f\n" % (haplotype, numerator, denominator))
                    
            # This heuristic is needed to account for edge cases in which a read
            # has a very large number of heterozygous sites, most of which match
            # a single phase. These can cause underflows in the denominator.
            try:
                denom = math.log10(denominator)
            except ValueError:
                denom = log10minfloat
            log_likelihood_ratios.append(numerator - denom)

        #sys.stderr.write("%s\n" % (log_likelihood_ratios))
        max_log_likelihood_ratio = max(log_likelihood_ratios)
        if max_log_likelihood_ratio == 0:
            phase = 'Unphased'
            max_phase = 0
        else:
            max_phase = log_likelihood_ratios.index(max_log_likelihood_ratio) + 1
#AB - this is a hard coded use of error rate instead of log likihood
#            if powlaw_modified(max_log_likelihood_ratio, a=powlaw_alpha, xmin=powlaw_xmin) <= error_rate_threshold:
                # if max_log_likelihood_ratio >= 10 ** (math.log10((error_rate_threshold / .5491)) / -2.45):
            phase = max_phase
#            else:
#                phase = 'Unphased'
        
        return [log_likelihood_ratios, phase, max_phase, max_log_likelihood_ratio]

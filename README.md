# HaplotagLR
## A tool for haplotagging long-read sequencing results.

HaplotagLR haplotags long reads based on existing, pre-phased haplotypes in VCF format.


## Dependencies
### All modes:
* HTSlib (https://www.htslib.org/)
* Python >= v3.7
* minimap2  (https://github.com/lh3/minimap2)
* numpy (https://numpy.org/)
* powerlaw (https://github.com/jeffalstott/powerlaw)
* pysam (https://github.com/pysam-developers/pysam)
* pyliftover (https://github.com/konstantint/pyliftover)
* requests (http://python-requests.org)
* samtools (https://github.com/samtools/samtools)

### Simulation mode
* pbsim2 (https://github.com/yukiteruono/pbsim2)


## Installation

We strongly recommend installing with conda, into a new environment:
```
conda create -n HaplotagLR_env -c conda-forge -c bioconda numpy pysam powerlaw pyliftover pbsim2 minimap2 requests samtools HaplotagLR python==3.7
```

Install with pip:
```
pip install HaplotagLR
```

Installation from the github repository is not recommended. However, if you must, follow the steps below:
1) git clone https://github.com/Boyle-Lab/HaplotagLR.git
2) cd HaplotagLR/
3) python3.7 -m pip install -e .


## Usage

```
HaplotagLR [-h] [--version] [-q] {haplotag,phasability,error_analysis} ...
```

HaplotagLR currently only offers haplotag mode, but may support more operations in future releases.


### Haplotag Mode
A tool for haplotagging individual long reads using pre-phased haplotypes

```
usage: HaplotagLR haplotag [-h] -v <VCF_FILE> -i <SAM/BAM/FASTQ>             
                           [-o </path/to/output>] [-r <REF_FASTA>]            
                           [-A <ASSEMBLY_NAME>] [-t <THREADS>] [-q] [-S] 
                           [-O {combined,phase_tagged,full}]
                           [-s <SAMPLE_NAME>] [-e EPSILON] [-c]
                           [-F FDR_THRESHOLD]                              
                           [--log_likelihood_threshold <MIN_LIKELIHOOD_RATIO>]
                           [--no_multcoeff]
```

#### Required Arguments
| Argument | Description |
|---|---|
| __-v <VCF_FILE>, --vcf <VCF_FILE>__ |Path to vcf file with haplotype information that will be used for haplotagging. (Must be in .vcf.gz format with tabix index in same folder. If .vcf file is provided, bgzip and tabix must be installed and available on PATH because HaplotagLR will attempt to convert it. EX: -v GM12878_haplotype.vcf.gz) |
| __-i <SAM/BAM/FASTQ>__ | Path to sequencing file (.fasta) or alignment file (.bam or .sam) of long reads that will be used for haplotagging. If either a .sam file is provided or an index is not found, .sam and .bam file will be sorted and indexed with SAMtools. Sorted.bam files should be in same directory as their index (.sorted.bam.bai). EX: -a data/minion_GM12878_run3.sorted.bam, -i minion_GM12878_run3.sam) Path to long read file in .fastq format that will be used for alignment and haplotagging (ex: -i minion_GM12878_run3.fastq). **** NOTE: the -r/--reference argument is REQUIRED if using input in fastq format! ****|

#### Optional Arguments
| Argument | Description |
|---|---|
| __-h, --help__ | Show help message and exit |
| __-o </path/to/output>, --output_directory_name </path/to/output_directory>__ | Output directory name. Name given to directory where results will be output. |
| __-r <REF_FASTA>, --reference <REF_FASTA>__ | Path to reference genome sequence file. REQUIRED if argument to -i a fastq file. |
| __-A <ASSEMBLY_NAME>, --reference_assembly <ASSEMBLY_NAME>__ | Assembly for the reference genome. EX: -A hg38. |
| __-t <THREADS>, --threads <THREADS>__ | Number of threads to use for mapping, sorting, and indexing steps. |
| __-q, --quiet__ | Output to stderr from subprocesses will be muted. |
| __-S, --silent__ | Output to stderr and stdout from subprocesses will be muted. |

#### Output Options
| Argument | Description |
|---|---|
| __-O {combined,phase_tagged,full}, --output_mode {combined,phase_tagged,full}__ | Specify whether/how phased, unphased, and nonphasable reads are printed to output. Modes available: combined: All reads will be written to a common output file. The phase tag (HP:i:N) can be used to extract maternal/paternal phased reads, unphased reads, and nonphasable reads. phase_tagged: Phased reads for both maternal and paternal phases will be written to a single output file, while unphased and nonphasable reads will be written to their own respective output files. full: Maternal, paternal, unphased, and nonphasable reads will be printed to separate output files. |
| __-s <SAMPLE_NAME>, --one_sample <SAMPLE_NAME>__ | Use the --one_sample option to haplotag a specific sample present in the input reads and vcf file. (-HG001) |

#### Statistical Options for Haplotag Mode
| Argument | Description |
|---|---|
| __-e EPSILON, --epsilon EPSILON__ | Use this value for the estimated sequencing error rate epsilon. This value will be used instead of calculating per-base error rates from quality scores, and in calculating the FDR threshold value rather than estimating the global mean sequencing error rate. Default = None. |
| __-c, --epsilon_from_quality_scores__ | Calculate per-base error rates directly from Phred scores in the BAM input. Default = False. |
| __-F FDR_THRESHOLD, --FDR_threshold FDR_THRESHOLD__ | Control the false discovery rate at the given value using a negative-binomial estimate of the number of haplotagging errors (N) given the average per-base sequencing error rate observed among all phaseable reads. Phased reads are sorted by their observed log-likelihood ratios and the bottom N*(1-FDR) reads will be reassigned to the "Unphased" set. Set this to zero to skip this step and return all haplotagging predictions. Default = 0.|
| __--log_likelihood_threshold <LOG_LIKELIHOOD_THRESHOLD>__ | Use a hard threshold on log-likelihood ratios when haplotagging reads. Results will only be printed for predicted haplotaggings with log-likelihood ratios equal to or greater than this threshold. Setting this to zero will cause all reads to be assigned to the phase to which they share the greatest number matches. Log-likelihood ratios will still be reported in the output in this case, but are not used for haplotagging decisions. |
| __--no_multcoeff__ | Do not apply the multinomial coefficient in the likelihood calculation. WARNING: The model will not output valid probabilities without this! Default=False (The multinomal coefficient will be used.) |

## Interpreting the Output
By default, HaplotagLR tags each record with several key:type:value tuples to encode the haplotagging decision and several values used in the tagging decision. These are written as BAM records to one or more output files, depending on invocation (See help for -O option). Specific tags added to each record are described below:

| Tag | Description |
|---|---|
| al | Alignment type, e.g., supplementary. Added if not already present in BAM record. |
| PS | Name of the overlapping phase set. |
| py | Ploidy number of overlapping phase set. |
| HS | Number of heterozygous variants overlapping read. |
| GP | Comma-delimited list of overlapping variants' position(s) relative to the read start. |
| PA | Phased alleles for all haplotypes overlapping read. Comma-delimited list of tuples. |
| OA | Observed allele(s) in sequenced read at heterozygous positions. May or may not match one of the values in PA! |
| PR | Comma-delimited list of Bayesian prior values. |
| LS | Comma-delimited list of Log-Likelihood-Ratios for each haplotype. |
| PC | Log-Likelihood-Ratio for assigned haplotag. |
| HP | Assigned haplotag. |

## Example Dataset
We provide a sample dataset and example usage [here](https://github.com/Boyle-Lab/HaplotagLR/tree/main/example_data)

## Citing HaplotagLR
The HaplotagLR algorithm and software release 1.0.3 are described in [pub link here](https://example.com). Please use the following citation if you use this software in your work:

HaplotagLR: an efficient algorithm for assigning haplotypic identity to long reads.
Monica J. Holmes, Babak Mahjour, Christopher Castro, Gregory A. Farnum, Adam G. Diehl, Alan P. Boyle.
2022. BioArxiv. [URL](http://example.com)

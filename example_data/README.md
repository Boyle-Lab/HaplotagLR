# LRphase_test_data
### Test dataset for LRphase

To run LRphase on the example dataset, first install LRphase according to the directions [here](https://github.com/Boyle-Lab/LRphase).

Next, download the example dataset (LRphase_Test_Data.tgz) to a suitable location from this link: [https://zenodo.org/record/7415834](https://zenodo.org/record/7415834) and unpack the tarball:
```
tar -xvzf LRphase_Test_Data.tgz
```

Working in the same location, issue the following command to produce combined phasing output for the example dataset:

```
LRphase phasing -v phased_variants.vcf.gz -i barcode01.Nanopore.sorted.bam -r hg38.no_alt.fa -A hg38.no_alt -s LIBD75
```

### Expected Results:

| File | MD5SUM |
| --- | --- |
| LIBD75.combined_phasing_results.bam | 35923f855f5b392f804e6f3b421e21f6 |

### Fastq data
Fastq data in barcode01.Nanopore.sorted.fastq.gz were generated from the bam file using bedtools bamtofastq. These can be used to generate equivalent output with the following command:

```
LRphase phasing -v phased_variants.vcf.gz -i barcode01.Nanopore.sorted.fastq.gz -r hg38.no_alt.fa -A hg38.no_alt -s LIBD75
```

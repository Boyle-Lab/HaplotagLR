from typing import Dict

hg19: Dict[str, str] = {}
hg38: Dict[str, str] = {}

hg38[
    'NA12877'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg38/small_variants/NA12877/NA12877.vcf.gz'
hg38[
    'NA12878'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz'

hg19[
    'NA12877'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg19/small_variants/NA12877/NA12877.vcf.gz'
hg19[
    'NA12878'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz'

hg38[
    'NA12891'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/10Xgenomics_ChromiumGenome_LongRanger2.1_09302016/NA12891_GRCh38/NA12891_GRCh38_phased_variants.vcf.gz'
hg38[
    'NA12892'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/10Xgenomics_ChromiumGenome_LongRanger2.1_09302016/NA12892_GRCh38/NA12892_GRCh38_phased_variants.vcf.gz'

hg19[
    'NA12891'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/10Xgenomics_ChromiumGenome_LongRanger2.1_09302016/NA12891_hg19/NA12891_hg19_phased_variants.vcf.gz'
hg19[
    'NA12892'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/10Xgenomics_ChromiumGenome_LongRanger2.1_09302016/NA12892_hg19/NA12892_hg19_phased_variants.vcf.gz'

hg38['reference_sequence_file_name'] = 'hg38.fa.gz'
hg19['reference_sequence_file_name'] = 'hg19.fa.gz'
hg38['reference_sequence'] = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
hg19['reference_sequence'] = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz'

hg38[
    'HG001'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
hg38[
    'HG002'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
hg38[
    'HG003'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
hg38[
    'HG004'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
hg38[
    'HG005'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz'
hg38[
    'HG006'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG006_NA24694_father/latest/GRCh38/HG006_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
hg38[
    'HG007'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG007_NA24695_mother/latest/GRCh38/HG007_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'

hg19[
    'HG001'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'
hg19[
    'HG002'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz'
hg19[
    'HG003'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz'
hg19[
    'HG004'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz'
hg19[
    'HG005'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/latest/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz'
hg19[
    'HG006'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG006_NA24694_father/latest/GRCh37/HG006_GIAB_GRCh37_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'
hg19[
    'HG007'] = 'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG007_NA24695_mother/latest/GRCh37/HG007_GIAB_GRCh37_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_highconf.vcf.gz'

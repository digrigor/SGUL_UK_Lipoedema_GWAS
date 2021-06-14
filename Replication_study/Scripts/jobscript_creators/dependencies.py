aggvcfs = "/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data"
funvcfs = "/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/functional_annotation/VEP"


exome_intervals = "/re_gecip/cardiovascular/postergaard/resources/Genome_reference_files/exome_calling_regions.v1.bed"
ref = "/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


gatk = "module load bio/GATK/4.1.4.1-GCCcore-8.3.0-Java-1.8;\n\ngatk"
bcftools = "module load bio/BCFtools/1.11-GCC-8.3.0;\n\nbcftools"
vcftools = "module load bio/VCFtools/0.1.16-foss-2018b-Perl-5.28.0;\n\nvcftools"
load_bcftools = "module load bio/BCFtools/1.11-GCC-8.3.0;\n\n"
load_vcftools = "module load bio/VCFtools/0.1.16-foss-2018b-Perl-5.28.0;\n\n"
bgzip = "module load bio/HTSlib/1.11-GCC-8.3.0;\n\nbgzip"
tabix = "module load bio/HTSlib/1.11-GCC-8.3.0;\n\ntabix"
rvtest = "module load bio/rvtests/2.1.0;\n\nrvtest"
vep = "module load lang/Perl/5.30.0-GCCcore-8.3.0; module load bio/VEP/99.1-foss-2019a-Perl-5.28.1; export PERL5LIB=$PERL5LIB:tools/apps/software/bio/VEP/99.1-foss-2019a-Perl-5.28.1/Plugins/; vep"
vep_filter = "module load bio/VEP/99.1-foss-2019a-Perl-5.28.1; filter_vep"
vcf2tsv = ". /resources/conda/miniconda3/etc/profile.d/conda.sh; conda activate py2_7_12pypi; python /re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/gel032020/vcf_melt_vcf2tsv.py"
plink2 = "module load bio/PLINK/2.00-devel-20200409-x86_64;\n\nplink2"
plink = "module load bio/PLINK/1.9b_4.1-x86_64;\n\n plink"




nthreads = 4

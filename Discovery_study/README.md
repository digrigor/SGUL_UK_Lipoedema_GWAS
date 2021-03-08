Discovery cohort: Genome-Wide Association Study
================
Dionysios Grigoriadis
2021-03-08

<style>
h1 {background: brown;color: white;padding-left: 7px;}
h2 {background: grey;color: white;padding-left: 7px;}
h3 {color: brown;}
</style>

FOR THE FULL VERSION OF THE ANALYSIS PLEASE OPEN THE HTML FILE IN THIS DIRECTORY

# Setup

We should first configure where Rstudio loads packages from, we just
need to run the setup script. We should also run the dependencies script
which contains paths of software or reference files we’re gonna use.

``` r
#Configure R studio packages location
library(ssh)
```

    ## Linking to libssh v0.8.6

``` r
run= F
stats3 = F

source('dependencies.R')

if(stats3==F){
  #session <- ssh_connect("dgrigori@stats3.sgul.ac.uk",)
}
```

    ## NULL

``` r
sr = paste0("~/Mimir",unlist(strsplit(getwd(),":"))[2])
```

Load Parameters + packages

``` r
library(rlang, quietly = T)
#install.packages("vctrs")
library(vctrs, quietly = T)
library(pillar, quietly = T)
library(magrittr, quietly = T)
library(dplyr, quietly = T)
library(knitr, quietly = T)
#library(stringi)
library(tidyselect, quietly = T)
library(plotly, quietly = T)
library(rprojroot, quietly = T)
library(rmarkdown, quietly = T)
library(DT, quietly = T)
#library(scales)
library(ggplot2, quietly = T)
library(kableExtra, quietly = T)
library(data.table, quietly = T)
#library(dplyr)
library(qqman, quietly = T)
library(CMplot, quietly = T)
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] CMplot_3.6.2      qqman_0.1.4       tidyselect_1.1.0  vctrs_0.3.6      
    ##  [5] rlang_0.4.10      ssh_0.7.0         data.table_1.13.6 kableExtra_1.3.1 
    ##  [9] scales_1.1.1      DT_0.17           rmarkdown_2.6     rprojroot_2.0.2  
    ## [13] plotly_4.9.3      ggplot2_3.3.3     stringi_1.5.3     knitr_1.31       
    ## [17] dplyr_1.0.4       magrittr_2.0.1    pillar_1.4.7     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.0.3    sys_3.4           tools_4.0.3       digest_0.6.27    
    ##  [5] packrat_0.5.0     jsonlite_1.7.2    evaluate_0.14     lifecycle_0.2.0  
    ##  [9] tibble_3.0.6      gtable_0.3.0      viridisLite_0.3.0 pkgconfig_2.0.3  
    ## [13] rstudioapi_0.13   yaml_2.2.1        xfun_0.20         stringr_1.4.0    
    ## [17] xml2_1.3.2        withr_2.4.1       httr_1.4.2        askpass_1.1      
    ## [21] generics_0.1.0    htmlwidgets_1.5.3 webshot_0.5.2     grid_4.0.3       
    ## [25] calibrate_1.7.7   glue_1.4.2        R6_2.5.0          purrr_0.3.4      
    ## [29] tidyr_1.1.2       MASS_7.3-53       credentials_1.3.0 ellipsis_0.3.1   
    ## [33] htmltools_0.5.1.1 rvest_0.3.6       colorspace_2.0-0  openssl_1.4.3    
    ## [37] lazyeval_0.2.2    munsell_0.5.0     crayon_1.4.0

# Lipoedema Batches Quality Controls

Peripheral blood from 148 lipoedema patients of white British ethnicity
were collected, DNA extracted and genotyped in two batches using
Illumina Infinium\_CoreExome-24\_v1-2 (Batch1) and
Infinium\_Core-24\_v1-2-a1 (Batch2) single nucleotide polymorphism (SNP)
chips, respectively.

## Files Management

``` r
sample_basename = "Input/Batch1_final/batch1_updated_all_prefixed"

dir = dirname(sample_basename)

qcdir = file.path(dir, "Qc/")

finalname = "lipoedema1_clean"

dir.create(dir)
dir.create(qcdir)
```

### Quality Control

Run PLINK Linux commands to:

**1)Generate the missing rate/calling rate, heterozygocity and sex
checking files **

./resources/software/plink1.9/plink –bfile
Input/Batch1\_final/batch1\_updated\_all\_prefixed –missing –het
–check-sex –out Input/Batch1\_final/Qc/batch1\_updated\_all\_prefixed

**2)Prune SNPs in relative LD so we can correctly calculate IBD **

./resources/software/plink1.9/plink –bfile
Input/Batch1\_final/batch1\_updated\_all\_prefixed –indep-pairwise 50 5
0.5 –out Input/Batch1\_final/Qc/batch1\_updated\_all\_prefixed

**3)IBD on pruned variants **

./resources/software/plink1.9/plink –bfile
Input/Batch1\_final/batch1\_updated\_all\_prefixed –maf 0.05 –genome
–out Input/Batch1\_final/Qc/batch1\_updated\_all\_prefixed –extract
Input/Batch1\_final/Qc/batch1\_updated\_all\_prefixed.prune.in

#### Input and Read data

Import all the files that plink generated into R. Then we can plot and
do some statistical analyses on them.

``` r
#Define the file names
imiss = file.path(qcdir,paste0(basename(sample_basename),'.imiss'))
het = file.path(qcdir,paste0(basename(sample_basename),'.het'))
sexcheck = file.path(qcdir,paste0(basename(sample_basename),'.sexcheck'))
ibd = file.path(qcdir,paste0(basename(sample_basename),'.genome'))

# Stop the script if the data file cannot be found in the current working directory:
stopifnot(file.exists(imiss))
stopifnot(file.exists(het))
stopifnot(file.exists(sexcheck))
stopifnot(file.exists(ibd))

#Load the actual data into R
miss_tab <- read.table(file=imiss, header=T, stringsAsFactors = FALSE)
het_tab <- read.table(file=het, header=T, stringsAsFactors = FALSE)
sexchk_tab <- read.table(file=sexcheck, header=T, stringsAsFactors = FALSE)
ibd_tab <- read.table(file=ibd, header=T, stringsAsFactors = FALSE)
```

#### Calculating Heterozygocity and Calling Rate per sample

Autosomal Heterozygocity = (Number of autosomal genotypes - Observed
number of homozygotes)/Number of autosomal genotypes  
Calling Rate = 1-(Number of missing genotype calls/ Number of valid
calls)  
Outliers: Samples with Heterozygocity \< mean-3SD heterozygocity calling
rate  
Samples with Heterozygocity \> mean+3SD heterozygocity calling rate  
Samples with Calling rate \< 0.97

``` r
#Calculate the heterozygocity rate
het_rate = (het_tab$N.NM. - het_tab$O.HOM.)/het_tab$N.NM.

#Calculate the genotype calling rate
cal_rate = 1-miss_tab$F_MISS

names(het_rate) = het_tab$IID
names(cal_rate) = miss_tab$IID

#Name the outliers
outliers = unique(c(names(cal_rate[cal_rate < 0.97]), names(het_rate[het_rate > mean(het_rate)+(3*sd(het_rate)) |                het_rate <  mean(het_rate)-(3*sd(het_rate))])))
outliers_rates = data.frame(Name=outliers, Call=cal_rate[outliers], Het=het_rate[outliers])

bline_outliers = unique(c(names(cal_rate[cal_rate < 0.975]), names(het_rate[het_rate > mean(het_rate)+(2*sd(het_rate))                  | het_rate < mean(het_rate)-(2*sd(het_rate))])))
bline_outliers_rates = data.frame(Name=bline_outliers, Call=cal_rate[bline_outliers], Het=het_rate[bline_outliers])

#Draw the heterozygocity to calling rate plot
hc_plot = ggplot(miss_tab, aes(x=het_rate, y=cal_rate)) + 
          geom_point(size=2, shape=21) +
          labs(x="Autosomal Heterozygocity Rate", y = "Calling Rate") + 
          theme_classic() +
          geom_vline(xintercept=mean(het_rate)+(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_vline(xintercept=mean(het_rate)-(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_hline(yintercept=0.97, colour='gray', linetype = "longdash", size=0.2) #+
          #geom_text(data=bline_outliers_rates, aes(x=Het, y=Call, label=Name), hjust=+1.2, vjust=0)
outliers = outliers

hc_plot
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Grey horizontal line: Calling Rate=0.97 threshold  
Green vertical line: ±3 SD deviation from the heterozygosity rate mean
of the samples threshold  
3 Outliers (will be excluded)

``` r
ibd_tab = ibd_tab[!(ibd_tab$IID1 %in% outliers | ibd_tab$IID2 %in% outliers),]
sexchk_tab = sexchk_tab[!sexchk_tab$IID %in% outliers,]
```

#### Check sex phenotypes from the genotype calls

Can our genotype calls correctly predict the sex of our samples?  
This is a QC step involving the calculation of the F inbreeding
coefficient of the X chromosome. Under sex-linkage, inbreeding only
causes an increase in the frequency of homozygous recessives in females.
There is no meaning to the inbreeding coefficient of a male since he is
hemizygous (i.e. he has only one X chromosome) and is, therefore, not
affected by inbreeding. Check
[here](http://www.genetic-genealogy.co.uk/supp/calc_inbreed.html) for
more information. Therefore, according to plink we are expecting:  
F \< 0.2 for females  
F \> 0.6 for males

``` r
#Get the F coefficient
sex_f = sexchk_tab$F
rownames(sexchk_tab)=paste0(sexchk_tab$FID,"-",sexchk_tab$IID)
#Create the plot
sex_f_color = ifelse(sexchk_tab$PEDSEX==sexchk_tab$SNPSEX, 'pink', 'grey')
barplot(sexchk_tab$F, col = sex_f_color, main='',
     xlab='F inbreeding coefficient', ylab='No. of samples')
legend('bottomright', legend=c('True Females', 'Females predicted as Males'), 
       col=c('pink','gray'), lwd=10)
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

According to this histogram, all of our samples had an F coefficient
less than 0.2 showing that our samples are reliably female samples.

Create a file of these IIDs and their FIDs which will be excluded. This
file will be used with the –remove plink command.

``` r
plink_exclude = miss_tab[miss_tab$IID %in% outliers,c('FID', 'IID')]
if(length(outliers)!=nrow(plink_exclude)){stop()}
#write.table(plink_exclude, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE)
```

### Create new plink clean set of files without the excluded samples:

We need to remove all four samples we decided to exclude from the
analysis:

./resources/software/plink1.9/plink –bfile
Input/Batch1\_final/batch1\_updated\_all\_prefixed –geno 0.05 –maf 0.01
–hwe 1e-6 midp –remove tempexcIDS.txt –make-bed –out
Input/Batch1\_final/lipoedema1\_clean; rm temp\*

## Lipoedema Batch 2

### Files management

``` r
sample_basename = "Input/Batch2_final/batch2_updated_all_prefixed"

dir = dirname(sample_basename)

qcdir = file.path(dir, "Qc/")

finalname = "lipoedema2_clean"

dir.create(dir)
```

    ## Warning in dir.create(dir): 'Input\Batch2_final' already exists

``` r
dir.create(qcdir)
```

    ## Warning in dir.create(qcdir): 'Input\Batch2_final\Qc' already exists

### Quality Control

Run PLINK Linux commands to:

**1)Generate the missing rate/calling rate, heterozygocity and sex
checking files **

./resources/software/plink1.9/plink –bfile
Input/Batch2\_final/batch2\_updated\_all\_prefixed –missing –het
–check-sex –out Input/Batch2\_final/Qc/batch2\_updated\_all\_prefixed

**2)Prune SNPs in relative LD so we can correctly calculate IBD **

./resources/software/plink1.9/plink –bfile
Input/Batch2\_final/batch2\_updated\_all\_prefixed –indep-pairwise 50 5
0.5 –out Input/Batch2\_final/Qc/batch2\_updated\_all\_prefixed

**3)IBD on pruned variants **

./resources/software/plink1.9/plink –bfile
Input/Batch2\_final/batch2\_updated\_all\_prefixed –maf 0.05 –genome
–out Input/Batch2\_final/Qc/batch2\_updated\_all\_prefixed –extract
Input/Batch2\_final/Qc/batch2\_updated\_all\_prefixed.prune.in

#### Input and Read data

Import all the files that plink generated into R. Then we can plot and
do some statistical analyses on them.

``` r
#Define the file names
imiss = file.path(qcdir,paste0(basename(sample_basename),'.imiss'))
het = file.path(qcdir,paste0(basename(sample_basename),'.het'))
sexcheck = file.path(qcdir,paste0(basename(sample_basename),'.sexcheck'))
ibd = file.path(qcdir,paste0(basename(sample_basename),'.genome'))

# Stop the script if the data file cannot be found in the current working directory:
stopifnot(file.exists(imiss))
stopifnot(file.exists(het))
stopifnot(file.exists(sexcheck))
stopifnot(file.exists(ibd))

#Load the actual data into R
miss_tab <- read.table(file=imiss, header=T, stringsAsFactors = FALSE)
het_tab <- read.table(file=het, header=T, stringsAsFactors = FALSE)
sexchk_tab <- read.table(file=sexcheck, header=T, stringsAsFactors = FALSE)
ibd_tab <- read.table(file=ibd, header=T, stringsAsFactors = FALSE)
```

#### Calculating Heterozygocity and Calling Rate per sample

Autosomal Heterozygocity = (Number of autosomal genotypes - Observed
number of homozygotes)/Number of autosomal genotypes  
Calling Rate = 1-(Number of missing genotype calls/ Number of valid
calls)  
Outliers: Samples with Heterozygocity \< mean-3SD heterozygocity calling
rate  
Samples with Heterozygocity \> mean+3SD heterozygocity calling rate  
Samples with Calling rate \< 0.97

``` r
#Calculate the heterozygocity rate
het_rate = (het_tab$N.NM. - het_tab$O.HOM.)/het_tab$N.NM.

#Calculate the genotype calling rate
cal_rate = 1-miss_tab$F_MISS

names(het_rate) = het_tab$IID
names(cal_rate) = miss_tab$IID

#Name the outliers
outliers = unique(c(names(cal_rate[cal_rate < 0.97]), names(het_rate[het_rate > mean(het_rate)+(3*sd(het_rate)) |                het_rate <  mean(het_rate)-(3*sd(het_rate))])))
outliers_rates = data.frame(Name=outliers, Call=cal_rate[outliers], Het=het_rate[outliers])

bline_outliers = unique(c(names(cal_rate[cal_rate < 0.975]), names(het_rate[het_rate > mean(het_rate)+(2*sd(het_rate))                  | het_rate < mean(het_rate)-(2*sd(het_rate))])))
bline_outliers_rates = data.frame(Name=bline_outliers, Call=cal_rate[bline_outliers], Het=het_rate[bline_outliers])

#Draw the heterozygocity to calling rate plot
hc_plot = ggplot(miss_tab, aes(x=het_rate, y=cal_rate)) + 
          geom_point(size=2, shape=21) +
          labs(x="Autosomal Heterozygocity Rate", y = "Calling Rate") + 
          theme_classic() +
          geom_vline(xintercept=mean(het_rate)+(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_vline(xintercept=mean(het_rate)-(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_hline(yintercept=0.97, colour='gray', linetype = "longdash", size=0.2) #+
          #geom_text(data=bline_outliers_rates, aes(x=Het, y=Call, label=Name), hjust=+1.2, vjust=0)
hc_plot
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
outliers = outliers
```

Grey horizontal line: Calling Rate=0.97 threshold Green vertical line:
±3 SD deviation from the heterozygosity rate mean of the samples
threshold 2 Outliers (will be excluded)

``` r
ibd_tab = ibd_tab[!(ibd_tab$IID1 %in% outliers | ibd_tab$IID2 %in% outliers),]
sexchk_tab = sexchk_tab[!sexchk_tab$IID %in% outliers,]
```

#### Check sex phenotypes from the genotype calls

Can our genotype calls correctly predict the sex of our samples?  
This is a QC step involving the calculation of the F inbreeding
coefficient of the X chromosome. Under sex-linkage, inbreeding only
causes an increase in the frequency of homozygous recessives in females.
There is no meaning to the inbreeding coefficient of a male since he is
hemizygous (i.e. he has only one X chromosome) and is, therefore, not
affected by inbreeding. Check
[here](http://www.genetic-genealogy.co.uk/supp/calc_inbreed.html) for
more information. Therefore, according to plink we are expecting:  
F \< 0.2 for females  
F \> 0.6 for males

``` r
#Get the F coefficient
sex_f = sexchk_tab$F
rownames(sexchk_tab)=paste0(sexchk_tab$FID,"-",sexchk_tab$IID)
#Create the plot
sex_f_color = ifelse(sexchk_tab$PEDSEX==sexchk_tab$SNPSEX, 'pink', 'grey')
barplot(sexchk_tab$F, col = sex_f_color, main='',
     xlab='F inbreeding coefficient', ylab='No. of samples')
legend('bottomright', legend=c('True Females', 'Females predicted as Males'), 
       col=c('pink','gray'), lwd=10)
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->
According to this histogram, all of our samples had an F coefficient
less than 0.2 showing that our samples are reliably female samples.

Create a file of these IIDs and their FIDs which will be excluded. This
file will be used with the –remove plink command.

``` r
plink_exclude = miss_tab[miss_tab$IID %in% outliers,c('FID', 'IID')]
if(length(outliers)!=nrow(plink_exclude)){stop()}
#write.table(plink_exclude, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE)
```

### Create new plink clean set of files without the excluded samples:

We need to remove all samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Input/Batch2\_final/batch2\_updated\_all\_prefixed –geno 0.05 –maf 0.01
–hwe 1e-6 midp –remove tempexcIDS.txt –make-bed –out
Input/Batch2\_final/lipoedema2\_clean; rm temp\*

# Merging Lipoedema Batches

### Files management

``` r
##Where are we going to run our analysis?
b1 = "Input/Batch1_final/lipoedema1_clean"
b2 = "Input/Batch2_final/lipoedema2_clean"

sample_basename = "Merged/Lipoedema_Batches_final/b1b2merged.no_at-cg.no_inc"

dir = dirname(sample_basename)
qcdir = file.path(dir, "Qc/")
qc2dir = file.path(qcdir, "Qc_arrays_inconsistency/")

finalname = "lipoedema_clean"

masterfile = "Input/Sample_info/kris_dx"
#lmaster = read.table(masterfile, sep="\t",stringsAsFactors = F)
#diagnoses <- lmaster$DIAGNOSIS.CONF.
#names(diagnoses) <- lmaster$D.Number

dir.create(dir, showWarnings = FALSE)
dir.create(qcdir, showWarnings = FALSE)
dir.create(qc2dir, showWarnings = FALSE)
```

### Load Input files

Load the bim and fam files of the plink filesets we want to merge

``` r
b1_bim = data.frame(fread(paste0(b1,".bim")), stringsAsFactors = FALSE)
b2_bim = data.frame(fread(paste0(b2,".bim")), stringsAsFactors = FALSE)

b1_fam = data.frame(fread(paste0(b1,".fam")), stringsAsFactors = FALSE)
b2_fam = data.frame(fread(paste0(b2,".fam")), stringsAsFactors = FALSE)

rownames(b1_bim)=b1_bim$V2
rownames(b2_bim)=b2_bim$V2

#Add a unique genomic position identifier (pid) with format [no of chromosome]chr[genomic coordinate]
b1_bim$pid = paste0(b1_bim$V1,'chr',b1_bim$V4)
b2_bim$pid = paste0(b2_bim$V1,'chr',b2_bim$V4)
```

### Keep only common variants

We need to only keep the common variants between the two batches.
Variants not covered in both filesets will generate unwanted biases and
batch effects.

**Keep in Batch1 only the variants which are shared between Batch2 and
Batch1**

less Input/Batch2\_final/lipoedema2\_clean.bim | awk ‘{print $2}’ \>
tempIDS.txt; ./resources/software/plink1.9/plink –bfile
Input/Batch1\_final/lipoedema1\_clean –extract tempIDS.txt –make-bed
–out Merged/Lipoedema\_Batches\_final/b1.shared

**Keep in Batch2 only the variants which are shared between Batch2 and
Batch1**

less Input/Batch1\_final/lipoedema1\_clean.bim | awk ‘{print $2}’ \>
tempIDS.txt; ./resources/software/plink1.9/plink –bfile
Input/Batch2\_final/lipoedema2\_clean –extract tempIDS.txt –make-bed
–out Merged/Lipoedema\_Batches\_final/b2.shared

## Merging

**1)Merge Batch2 and Filtered Batch1 (This step will fail due to strand
flipping issues but it will create the usefull missnp file which will be
used for strand flipping.)**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b2.shared –bmerge
Merged/Lipoedema\_Batches\_final/b1.shared

**2)Use the missnp file to filter flip the problematic snps in batch1
and then redo the merging.**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1.shared –flip plink.missnp –make-bed
–out Merged/Lipoedema\_Batches\_final/b1.shared.flipped

Now we can actually merge batch2 with the flipped version of batch1
**3)Merge the flipped batch1 with the batch2 data-set **  
In this step **WE ALSO REMOVE NON ACGT SNPS from the merged dataset**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b2.shared –bmerge
Merged/Lipoedema\_Batches\_final/b1.shared.flipped –merge-mode 1
–snps-only just-acgt –make-bed –out
Merged/Lipoedema\_Batches\_final/b1b2merged

Merge is now completed.

## Quality Control

### T -\> A, A -\> T, G -\> C, C -\> G SNPs removal

T -\> A, A -\> T, G -\> C, C -\> G SNPs are particularly hard to resolve
when merging datasets. In this case we are going to remove them from the
merged dataset.

    ## 652 T -> A, A -> T, G -> C, C -> G SNPs are being removed:

    ## kgp776941, rs4646515, rs307354, rs2843130, rs11121382, rs848578, kgp468659, GA000342, rs665691, rs11247703, rs501832, rs837398, rs2297812, rs10788882, rs671108, rs7523104, kgp12171305, rs6588109, rs549148, rs10489626, rs17131479, kgp7261306, rs12730292, rs3738573, rs12144715, rs2503243, rs12749687, rs3843306, rs12030843, kgp7546508, rs1520668, rs2786529, rs309087, rs1124427, rs521734, kgp3155172, kgp8762234, kgp6470243, rs2488429, kgp1092768, kgp10847070, rs3136701, kgp1727928, kgp2698943, rs1212352, rs6666258, rs11537583, rs4147595, rs2236866, rs4987389, kgp6541618, kgp1999495, kgp7005320, rs17649050, rs859715, rs609521, rs12568050, rs832174, rs6693954, rs2243158, rs12130212, rs7524203, rs4149230, rs10916025, rs369909, rs514230, rs10926556, rs3102947, rs478222, rs4666002, kgp522146, kgp1097343, rs1031261, kgp3277363, rs891706, rs3816182, rs17035884, rs12620679, rs4233949, rs12713268, rs10184151, rs2234500, rs2916519, rs867529, rs2314398, rs3748930, rs13151, rs1402467, rs2276561, rs11123170, rs4513299, rs10207628, rs13025959, rs1050900, kgp5104289, rs6735208, rs6736104, rs6741949, rs1521527, rs1993308, kgp11638042, rs11757, rs2037799, rs12476147, rs10184895, rs7582694, rs2176528, rs12693898, rs17448457, rs4675644, rs2229571, rs12987009, rs1983210, rs2047136, rs16824283, rs13395911, kgp7355471, rs1042640, kgp7359081, rs4305276, rs1003972, rs402675, rs2279975, rs271066, rs2290159, rs7640928, rs4858649, rs7630967, kgp6063557, rs7372943, rs2271087, rs1468542, rs1005678, rs10490778, rs7433808, rs13069000, rs13069049, rs11713321, rs12486865, rs1479371, rs2276872, rs9811423, kgp11508564, rs4677948, rs10512627, kgp12182453, rs8649, rs5013525, rs4894410, rs9832727, kgp2007953, rs4149496, rs10935807, rs6798100, rs7628548, rs11920090, rs7430671, rs10937329, rs17505102, rs6815464, rs11734132, rs3775948, rs11733056, rs1972463, rs13101978, rs10014198, rs279845, rs726967, rs4148271, kgp6788842, rs3775782, rs13133166, rs344141, rs931605, rs578017, rs10222732, rs10019009, rs2725220, kgp3414215, kgp2808767, rs17219837, rs4147541, 401070, kgp3087994, kgp11897209, rs2306986, kgp1506779, kgp1489278, rs2132845, rs10028213, rs1466662, rs17046216, rs7697970, rs4241816, rs1593, rs1062547, rs7447815, rs4147775, rs7702187, rs3822410, rs4701523, rs12523160, rs1063499, rs7709645, rs2305962, rs164572, rs17663555, rs961098, kgp3636818, kgp1037316, rs17513503, rs7705033, kgp775040, rs274547, rs2548966, rs744674, rs12653537, rs3792796, rs13358864, rs2278255, rs2070998, rs1993554, rs3798713, kgp641225, rs7747752, rs2010190, rs198843, rs9467750, kgp5879917, rs454182, rs7743296, rs209151, rs3094548, rs1010408, rs3094736, rs1611737, rs9261578, rs1264540, rs2534823, rs3132613, rs4713422, rs2233967, rs130078, rs34282959, rs9266399, rs7775117, rs2523453, rs2523685, rs3828890, rs3853601, rs9332739, rs6455, rs9267845, rs1041885, rs35265698, rs9271300, kgp26489400, kgp12191290, kgp1963903, kgp7884518, kgp2881793, rs4713605, rs6913896, rs734181, rs9394159, rs831510, kgp5152142, rs4714750, kgp597715, rs2206271, rs9474433, rs561930, rs9446917, rs13202860, rs11755527, rs211221, kgp9903653, rs7760535, rs619203, rs6569474, rs9478858, kgp10632871, exm592133, rs662138, rs3818678, kgp7801610, rs3127594, rs2504916, kgp648939, rs12214416, rs7772437, rs4252105, GA006079, rs10275044, rs10277115, rs10245377, rs852499, rs6967241, kgp3981052, kgp4690887, kgp4722850, rs324981, rs1065647, rs28637820, rs884103, rs10245531, rs4148830, kgp5355388, rs2888611, rs17064, kgp1123831, rs705382, rs10487131, kgp3815960, rs7801803, rs4727338, kgp12362128, rs4729656, rs41521, rs6958498, rs1404697, rs4148686, rs4609139, rs2598291, rs194150, rs2240395, rs7781853, rs6944935, rs10282008, rs10263087, rs3800855, rs2407314, rs930557, rs6984305, rs1484642, rs12458, rs2460338, rs437548, rs7833268, rs2614091, rs2439302, rs12677427, rs2555588, rs2380571, rs10100101, rs2167801, rs7007970, rs2875931, rs2001945, kgp1315012, rs9644251, rs6558295, rs6474795, rs3731211, rs4978053, rs7468614, rs3043, kgp3137687, rs914715, kgp10939539, rs10780660, rs10868406, rs10116326, rs10118939, rs10819937, rs9409038, rs3824479, rs1929842, kgp2407207, kgp27322098, rs3204145, rs1327796, rs3818764, rs7852872, rs7039505, rs7851693, rs2073873, rs11103473, rs12722605, rs2501677, rs4748140, rs4748309, rs12258967, rs4747442, rs1148186, rs4590817, kgp5504085, rs1871452, rs754466, kgp12472133, kgp10638206, rs7923151, rs3758383, kgp10852535, kgp7668749, rs1853205, rs28399513, kgp537287, rs17222723, kgp7613533, rs521674, rs2065779, rs12413624, rs1999628, rs7923707, kgp5910402, rs4752781, rs689, rs163182, rs1026231, rs10500633, rs2722769, rs11022131, rs11022257, rs11023332, rs1793003, rs391317, rs2074040, rs4148904, rs12285276, rs2070852, rs10501320, rs1945213, rs174479, rs2276299, kgp5311770, rs504915, rs10896027, kgp4375382, rs2276288, rs597480, rs7396702, rs10830962, rs10830989, rs11021436, rs633185, rs516091, rs4938016, rs17118846, rs12225230, rs4938642, rs11218350, rs11222084, rs4765905, rs4765913, rs11053646, rs4149192, rs4149197, rs10841684, kgp277759, rs5488, kgp6443294, kgp11801419, rs864360, rs10842262, rs12819160, rs7953528, rs3759302, rs11564148, rs736825, rs2638315, rs941207, rs1495377, rs11178997, rs10879392, rs2365919, rs7300275, rs12423247, rs4584635, rs4766578, rs10849893, rs2388082, rs12323063, rs3812883, rs3736830, rs1063181, rs287398, rs9543067, rs1751051, rs279942, kgp8633405, rs12873360, rs4773279, kgp9114043, rs2273394, rs12879346, rs9323124, rs12436916, rs1061108, rs1997896, rs2025009, rs2239557, rs2232700, rs61280460, rs2069590, rs11628318, rs3178152, rs2273703, rs4906695, rs2178004, rs11635685, rs7166580, rs9325, kgp9514200, rs3204689, rs4775041, rs8034331, rs3743171, rs4776970, rs17279860, rs1874953, rs17875563, rs2521501, rs7165226, GA030783, kgp3881218, rs7168353, rs970843, kgp11779065, rs1005190, rs323043, rs11076863, rs11643447, rs4148344, kgp799675, rs212090, kgp16452123, kgp5026912, rs7190447, rs9924528, rs11649653, kgp4393797, rs12447652, rs2917666, rs4788815, rs16973585, rs7186908, kgp9244206, kgp961638, kgp6232103, kgp6429313, rs11149631, rs1060253, rs1056707, rs12051548, rs12449580, rs4151120, rs4462665, rs4985726, kgp7846805, rs2228100, rs9908604, rs8073060, rs11653310, rs16971217, rs2293152, rs17651507, rs199515, rs12452766, rs11868380, rs12451254, rs11658329, rs11653989, rs7501742, rs2853533, rs8090011, kgp9021777, rs1941059, kgp1984670, rs567458, rs10853545, rs125555, rs1037757, rs8259, rs2074442, rs1026466, rs2446195, rs428453, rs2230199, rs2161468, rs10419865, rs2738446, rs714772, rs4646529, rs2239367, rs10418069, kgp6160768, rs12977516, rs2305802, rs773930, rs401502, rs11671559, rs2967456, rs28493229, kgp5448752, rs11547806, rs365745, rs519113, rs7259004, rs6966, rs1064202, rs296365, rs2307279, rs10417472, rs3829655, rs8105668, rs4813043, rs1053783, rs6040399, rs2250889, rs12617, kgp3493774, rs2869343, rs7271624, rs6026990, rs2032297, rs198887, rs9983044, rs25678, rs2839398, kgp10445119, rs4148111, rs6586345, kgp11497610, kgp6830733, kgp2318663, rs17004785, kgp23307465, rs165656, kgp11396666, kgp1866281, rs4820255, rs228941, rs229526, rs6002455, kgp3971415, rs2071885, rs738409, rs6007344, kgp8465547, kgp9006041, kgp6121392, kgp2049186, rs3087965, rs17312163, rs3810682, rs10521351, rs1017682, rs1494481, rs1180895, rs1734791, kgp22775108

Use plink to remove these SNPs

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged –exclude tempIDS.txt
–make-bed –out Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg

### Remove SNPs that are incosistent between the same samples genotyped in different batches

**1) First we extract only good quality SNPs from the merged dataset
(maf, geno, hwe filters) and export the resulting dataset as vcf file**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg –maf 0.01 –geno
0.05 –hwe 1e-6 midp –recode vcf –out
Merged/Lipoedema\_Batches\_final/Qc/Qc\_arrays\_inconsistency/b1b2merged.no\_at-cg.maf.geno.hwe

**2) Actual check for inconsistencies between the 2 lipoedema arrays**

To avoid batch effect generated by genotyping the lipoedema samples in
two slightly different SNP arrays, 22 samples were genotyped in both
batches.

``` r
#Load the VCF
vcf = data.frame(fread(file.path(qc2dir,paste0("b1b2merged.no_at-cg.maf.geno.hwe",".vcf"))),stringsAsFactors = F)

b1sams = colnames(vcf)[startsWith(colnames(vcf),"B1_")]
b2sams = colnames(vcf)[startsWith(colnames(vcf),"B2_")]

names(b1sams) = unlist(lapply(colnames(vcf)[startsWith(colnames(vcf),"B1_")], function(x) unlist(strsplit(x, "_"))[3]))
names(b2sams) = unlist(lapply(colnames(vcf)[startsWith(colnames(vcf),"B2_")], function(x) unlist(strsplit(x, "_"))[3]))

b1_shared = c(b1sams[which(names(b1sams) %in% names(b2sams))],b1sams["D116479"])
b2_shared = c(b2sams[which(names(b2sams) %in% names(b1sams))],b2sams["D116497"])

#Check for inconsistent IDs per samples
incids = unlist(lapply(c(1:length(b1_shared)), function(x) vcf[(vcf[,b1_shared[x]]!=vcf[,b2_shared[x]] & vcf[,b1_shared[x]]!="./." & vcf[,b2_shared[x]]!="./."),]$ID))

incsnps = unique(unlist(incids))

#Write all the inconsistent SNPs on a txt file which we will use to exclude them
#write.table(incsnps, paste0(qcdir,'/incIDS.txt'), quote=FALSE, row.names = FALSE, col.names = FALSE)
#write.table(incsnps, 'tempIDS.txt', quote=FALSE, row.names = FALSE, col.names = FALSE)

#Write inconsistent tables for illumina
verincsnps = names(which(sort(table(unlist(incids)))>3))
vcfb1 = vcf[vcf$ID %in% verincsnps, c(colnames(vcf)[1:9],b1_shared)]
colnames(vcfb1) = unlist(lapply(colnames(vcfb1), function(x) tail(unlist(strsplit(x,"_")),n=1)))
vcfb2 = vcf[vcf$ID %in% verincsnps, c(colnames(vcf)[1:9],b2_shared)]
colnames(vcfb2) = unlist(lapply(colnames(vcfb2), function(x) tail(unlist(strsplit(x,"_")),n=1)))
#write.csv(vcfb1, file.path(qc2dir,"inconsistent_snps_b1.csv"), quote = F, row.names = F)
#write.csv(vcfb2, file.path(qc2dir,"inconsistent_snps_b2.csv"), quote = F, row.names = F)

cat(paste0(length(unique(c(vcfb1$ID,vcfb2$ID)))," SNPs showing inconsistency between the 2 arrays were removed: "))
```

    ## 4 SNPs showing inconsistency between the 2 arrays were removed:

``` r
cat(paste(unique(c(vcfb1$ID,vcfb2$ID)),collapse=", "))
```

    ## rs2382444, rs10904494, rs8015827, rs2248797

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg –exclude
tempIDS.txt –make-bed –out
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc

### Exclude Badly phenotyped lipoedema samples from the analysis

``` r
lipo_to_keep = read.table("Input/Sample_info/lipoedema_samples_touse.txt",stringsAsFactors = F, header=T)

lipo_excludes = mfam[!mfam$V2 %in% lipo_to_keep$DNumber,c(1,2)]

#write.table(lipo_excludes, file.path(dir,'lipo_excludes.txt'), sep='\t', row.names = FALSE, quote=FALSE)

new_sample_basename = "Merged/Lipoedema_Batches_final/b1b2merged.no_at-cg.no_inc.true_lipo"
```

We need to remove all samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc –geno 0.05
–maf 0.01 –hwe 1e-6 midp –remove
Merged/Lipoedema\_Batches\_final/lipo\_excludes.txt –make-bed –out
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc.true\_lipo;
rm temp\*

### Set a new basename

``` r
sample_basename = new_sample_basename
```

### SNP and Samples Quality Control

**1)Generate the missing rate/calling rate, heterozygocity and sex
checking files **

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc.true\_lipo
–missing –het –check-sex –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo

**2)Prune SNPs in relative LD so we can correctly calculate IBD **

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc.true\_lipo
–indep-pairwise 50 5 0.5 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo

**3)IBD on pruned variants **

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc.true\_lipo
–maf 0.01 –genome –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo
–extract
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.prune.in

#### Input and Read data

Import all the files that plink generated into R. Then we can plot and
do some statistical analyses on them.

``` r
#Define the file names
imiss = file.path(qcdir,paste0(basename(sample_basename),'.imiss'))
het = file.path(qcdir,paste0(basename(sample_basename),'.het'))
sexcheck = file.path(qcdir,paste0(basename(sample_basename),'.sexcheck'))
ibd = file.path(qcdir,paste0(basename(sample_basename),'.genome'))

# Stop the script if the data file cannot be found in the current working directory:
stopifnot(file.exists(imiss))
stopifnot(file.exists(het))
stopifnot(file.exists(sexcheck))
stopifnot(file.exists(ibd))

#Load the actual data into R
miss_tab <- read.table(file=imiss, header=T, stringsAsFactors = FALSE)
het_tab <- read.table(file=het, header=T, stringsAsFactors = FALSE)
sexchk_tab <- read.table(file=sexcheck, header=T, stringsAsFactors = FALSE)
ibd_tab <- read.table(file=ibd, header=T, stringsAsFactors = FALSE)

ibd_tab$ID1 = paste(ibd_tab$FID1,ibd_tab$IID1,sep="-")
ibd_tab$ID2 = paste(ibd_tab$FID2,ibd_tab$IID2,sep="-")
sexchk_tab$ID = paste(sexchk_tab$FID,sexchk_tab$IID,sep="-")
```

# QC Measurements

## Calculating Heterozygocity and Calling Rate per sample

Autosomal Heterozygocity = (Number of autosomal genotypes - Observed
number of homozygotes)/Number of autosomal genotypes  
Calling Rate = 1-(Number of missing genotype calls/ Number of valid
calls)  
Outliers: Samples with Heterozygocity \< mean-3SD heterozygocity calling
rate  
Samples with Heterozygocity \> mean+3SD heterozygocity calling rate  
Samples with Calling rate \< 0.97

``` r
#Calculate the heterozygocity rate
het_rate = (het_tab$N.NM. - het_tab$O.HOM.)/het_tab$N.NM.

#Calculate the genotype calling rate
cal_rate = 1-miss_tab$F_MISS

names(het_rate) = paste(het_tab$FID,het_tab$IID,sep="-")
names(cal_rate) = paste(miss_tab$FID,miss_tab$IID,sep="-")

ibd_tab$CAL1 = cal_rate[ibd_tab$ID1]
ibd_tab$CAL2 = cal_rate[ibd_tab$ID2]


#Name the outliers
outliers = unique(c(names(cal_rate[cal_rate < 0.97]), names(het_rate[het_rate > mean(het_rate)+(3*sd(het_rate)) |                het_rate <  mean(het_rate)-(3*sd(het_rate))])))
outliers_rates = data.frame(Name=outliers, Call=cal_rate[outliers], Het=het_rate[outliers])

bline_outliers = unique(c(names(cal_rate[cal_rate < 0.975]), names(het_rate[het_rate > mean(het_rate)+(2*sd(het_rate))                  | het_rate < mean(het_rate)-(2*sd(het_rate))])))
bline_outliers_rates = data.frame(Name=bline_outliers, Call=cal_rate[bline_outliers], Het=het_rate[bline_outliers])

#Draw the heterozygocity to calling rate plot
hc_plot = ggplot(miss_tab, aes(x=het_rate, y=cal_rate)) + 
          geom_point(size=2, shape=21) +
          labs(x="Autosomal Heterozygocity Rate", y = "Calling Rate") + 
          theme_classic() +
          geom_vline(xintercept=mean(het_rate)+(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_vline(xintercept=mean(het_rate)-(3*sd(het_rate)), colour='green', linetype = "longdash", size=0.2) +
          geom_hline(yintercept=0.97, colour='gray', linetype = "longdash", size=0.2) #+
          #geom_text(data=outliers_rates, aes(x=Het, y=Call, label=Name), hjust=0, vjust=+1.5, size=3)
hc_plot
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-86-1.png)<!-- -->

``` r
outliers = outliers
```

Grey horizontal line: Calling Rate=0.97 threshold Green vertical line:
±3 SD deviation from the heterozygosity rate mean of the samples
threshold 2 Outliers (will be excluded)

``` r
ibd_tab = ibd_tab[!(ibd_tab$ID1 %in% outliers | ibd_tab$ID2 %in% outliers),]
sexchk_tab = sexchk_tab[!sexchk_tab$ID %in% outliers,]
```

## Check sex phenotypes from the genotype calls

Can our genotype calls correctly predict the sex of our samples?  
This is a QC step involving the calculation of the F inbreeding
coefficient of the X chromosome. Under sex-linkage, inbreeding only
causes an increase in the frequency of homozygous recessives in females.
There is no meaning to the inbreeding coefficient of a male since he is
hemizygous (i.e. he has only one X chromosome) and is, therefore, not
affected by inbreeding. Check
[here](http://www.genetic-genealogy.co.uk/supp/calc_inbreed.html) for
more information. Therefore, according to plink we are expecting:  
F \< 0.2 for females  
F \> 0.6 for males

``` r
#Get the F coefficient
sex_f = sexchk_tab$F
rownames(sexchk_tab)=paste0(sexchk_tab$FID,"-",sexchk_tab$IID)
#Create the plot
sex_f_color = ifelse(sexchk_tab$PEDSEX==sexchk_tab$SNPSEX, 'pink', 'grey')
barplot(sexchk_tab$F, col = sex_f_color, main='',
     xlab='F inbreeding coefficient', ylab='No. of samples')
legend('bottomright', legend=c('True Females', 'Females predicted as Males'), 
       col=c('pink','gray'), lwd=10)
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-89-1.png)<!-- -->
According to this histogram, all of our samples had an F coefficient
less than 0.2 showing that our samples are reliably female samples.

## Check relatedness from the genotype calls by measuring Identity-by-descent (IBD)

Here we plot the IBD of each sample pair in our samples (Excluding the
previously pointed outliers).  
In our data-set there are some related samples which we need to exclude.
We therefore expect an IBD = 0 for most of the samples but not for these
related ones.

``` r
#Plot the IBD of all samples without the previously excluded ones
barplot(ibd_tab$PI_HAT, main="IBD plot before relateds exclusion", xlab = 'Samples', ylab='Identity-by-descent (IBD)', ylim=c(0,1))
abline(h = 0.05, lty = 2, col='red')
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-90-1.png)<!-- -->

``` r
ibd_tab$CAL1 = cal_rate[ibd_tab$ID1]
ibd_tab$CAL2 = cal_rate[ibd_tab$ID2]
ibd_tab$ID1 = paste0(ibd_tab$FID1,"-",ibd_tab$IID1)
ibd_tab$ID2 = paste0(ibd_tab$FID2,"-",ibd_tab$IID2)

#write.table(ibd_tab,file.path(qcdir,paste0(basename(sample_basename),'_before.genome')), quote = F, row.names = F)
#write.table(subset(ibd_tab,ibd_tab$PI_HAT>0.05),file.path(qcdir,"ibdcheck.txt"), quote = F, row.names = F)
```

Relatedness between all sample pairs in the cohort was inferred by
calculating identity by descent. In sample-pairs with PI\_HAT\>0.05, the
sample with less profound lipoedema characteristics (for cases) and/or
lower genotyping calling rate was excluded:

``` r
toexc = read.table(paste0(qcdir,"/","ibdcheck_to_exclude.txt"),stringsAsFactors = F,header=F)
colnames(toexc) = c("ID")

ibd_tab = ibd_tab[!(ibd_tab$ID1 %in% toexc$ID | ibd_tab$ID2 %in% toexc$ID),]
barplot(ibd_tab$PI_HAT, main="IBD plot AFTER relateds exclusion", xlab = 'Samples', ylab='Identity-by-descent (IBD)', ylim=c(0,1))
abline(h = 0.05, lty = 2, col='red')
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-91-1.png)<!-- -->

``` r
ibd_exc=toexc$ID
```

    ## 27 outliers will be excluded due to relatedness

``` r
#In the previously defined outlier list, add the related samples
outliers = c(outliers, unique(ibd_exc))

toexcdf_outs = data.frame(FID = unlist(lapply(outliers, function(x) unlist(strsplit(x, "-"))[1])), IID = unlist(lapply(outliers, function(x) unlist(strsplit(x, "-"))[2])))

toexcdf = toexcdf_outs
```

Perfect\!

Create a file of these IIDs and their FIDs which will be excluded. This
file will be used with the –remove plink command.

``` r
#write.table(toexcdf, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE)
```

## PCA Analysis

### Create new plink clean set of files without the excluded samples:

We need to remove all samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged.no\_at-cg.no\_inc.true\_lipo
–remove tempexcIDS.txt –make-bed –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd;
rm temp*; rm toexc*;

### LD prune and remove any missing calls from the new clean data-set (Like we did in the beginning of the analysis)

***1)Prune SNPs in relative LD and remove missing calls so we can
correctly calculate Principal Components later ***

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd
–indep-pairwise 50 5 0.5 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldprune

**2) Remove the pruned variants from the dataset **

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd
–maf 0.01 –geno 0.05 –genome –extract
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldprune.prune.in
–make-bed -out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldpruned

### Calculation of principal components (PCA Plot) of our clean and pruned Data-Set

We will generate the PCA plot to check how the samples are clustered
based on their genotype data. Any batch effects such as ethnicity will
be revealed. We’ll be using GCTA software for this (check Dependencies.R
script).

First, we need to convert our plink files to GCTA grm format

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldpruned
–make-grm –thread-num 10 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldpruned

Now we can use the grm files we created to perform the GCTA PC analysis.
The output eigenvec and eigenval files are the ones we need for the PCA
plot

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldpruned
–pca 4 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd.ldpruned

\#\#\#PCA Plot on the LD pruned data-set

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)
pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)
pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))
colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
pcs$group = NA
pcs[startsWith(pcs$FID,"B1_"),]$group = "Batch1"
pcs[startsWith(pcs$FID,"B2_"),]$group = "Batch2"

#PC1 vs PC2
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC2") +
          theme_classic()

#PC2 vs PC3
pca2 = ggplot(pcs, aes(x=PC2, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC2", y = "PC3") +
          theme_classic()

#PC1 vs PC3
pca3 = ggplot(pcs, aes(x=PC1, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC3") +
          theme_classic()

#Table with EigenValues and all the PCA plots
pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues') %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Eigenvalues

</caption>

<thead>

<tr>

<th style="text-align:right;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1.26552

</td>

</tr>

<tr>

<td style="text-align:right;">

1.13110

</td>

</tr>

<tr>

<td style="text-align:right;">

1.09867

</td>

</tr>

<tr>

<td style="text-align:right;">

1.09413

</td>

</tr>

<tr>

<td style="text-align:right;">

1.09061

</td>

</tr>

</tbody>

</table>

``` r
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-110-1.png)<!-- -->

``` r
pca2
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-110-2.png)<!-- -->

``` r
pca3
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-110-3.png)<!-- -->

There are outliers in our PCA plots which might need to be excluded. To
make more sense of it we’re going to merge our dataset with the HapMap
dataset and create a new PCA which will show any ancestry outliers and
their origin.

### Merge our data with three Hapmap populations; CEU; Yoruban and CHB.

**PCA analysis of GWAS & HapMap samples** We expect all our samples to
cluster together with the Central European samples.

**1)Extract HapMap SNPs **

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd
–extract ./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.no-at-cg-snps.txt
–make-bed –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps

**2)We will merge the HapMap and Batch2 datasets but multi-allelic SNPs
will be created: If an allele is the major in one dataset and the minor
in the other dataset, then plink marks them as multi-allelic SNPs which
we will remove.**

**Initially merge the HapMap and Batch2 datasets:**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed

**Remove the multi-allelic SNPs**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps
–exclude plink.missnp –make-bed –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed
–out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Prune SNPs in relative LD and remove missing calls so we can correctly
calculate Principal Components later**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA
–indep-pairwise 50 5 0.5 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Remove the pruned variants from the dataset**

./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA
–maf 0.01 –geno 0.05 –extract
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.prune.in
–make-bed -out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Convert our plink files to GCTA grm format **

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–make-grm –autosome –thread-num 10 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Now we can use the grm files we created to perform the GCTA PC
analysis. The output eigenvec and eigenval files are the ones we need
for the PCA plot **

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–pca 4 –out
Merged/Lipoedema\_Batches\_final/Qc/b1b2merged.no\_at-cg.no\_inc.true\_lipo.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)
pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)
colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
pcs$group = NA
pcs[startsWith(pcs$FID,"B1_"),]$group = "Lipoedema_Batch1"
pcs[startsWith(pcs$FID,"B2_"),]$group = "Lipoedema_Batch2"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"1"),]$group = "CEU"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"NA"),]$group = "CHB"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"Y"),]$group = 'YRI'
pcs$group = factor(pcs$group)


pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))


#pc_vals = read.table(paste0(plink_dir,out_suffix,'_hapmap-snps2_for_PCA_pruned.eigenval'), header=FALSE, stringsAsFactors = FALSE)

#colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
#colnames(pc_vals)=c('EigenValues')

#PC1 vs PC2
#pcs = pcs[pcs$group %in% c("Lipoedema_Batch1","Lipoedema_Batch2"),]
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group, label=paste0(FID,"-",IID))) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC2") + geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          theme_classic()
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-135-1.png)<!-- -->

``` r
pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues')  %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Eigenvalues

</caption>

<thead>

<tr>

<th style="text-align:right;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

41.85960

</td>

</tr>

<tr>

<td style="text-align:right;">

35.75280

</td>

</tr>

<tr>

<td style="text-align:right;">

1.86259

</td>

</tr>

<tr>

<td style="text-align:right;">

1.51920

</td>

</tr>

<tr>

<td style="text-align:right;">

1.36796

</td>

</tr>

</tbody>

</table>

``` r
#gg = ggplotly(pca1)
#gg
```

    ## 5  PCA Outliers will be excluded.

These are outliers around the European population and therefeore we
excluded them.

## Create the final clean set of files:

Save the ids that need to be removed:

We need to remove all the samples and variants we decided to exclude
from the analysis:

cat Merged/Lipoedema\_Batches\_final/Qc/incIDS.txt \>
Merged/Lipoedema\_Batches\_final/Qc/incatcgIDS.txt; cat
Merged/Lipoedema\_Batches\_final/Qc/atcgIDS.txt \>\>
Merged/Lipoedema\_Batches\_final/Qc/incatcgIDS.txt;
./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/b1b2merged –exclude
Merged/Lipoedema\_Batches\_final/Qc/incatcgIDS.txt –remove
Merged/Lipoedema\_Batches\_final/Qc/removeidsfinal.txt –maf 0.01 –geno
0.05 –hwe 1e-6 midp –make-bed –out
Merged/Lipoedema\_Batches\_final/lipoedema\_clean

# Controls - Understanding Society Data

5849 female samples enrolled in the Understanding Society UK study
(Understanding Society UK) and genotyped using HumanCoreExome-12\_v1.0
were used as controls (European Genome-phenome Archive ID:
EGAD00010000890).

\#Here we need to specify some locations for the analysis

``` r
sample_basename = "Input/USociety_final/USociety"

dir = dirname(sample_basename)

qcdir = file.path(dir, "Qc/")

finalname = "USociety_clean"

dir.create(dir)
```

    ## Warning in dir.create(dir): 'Input\USociety_final' already exists

``` r
dir.create(qcdir)
```

    ## Warning in dir.create(qcdir): 'Input\USociety_final\Qc' already exists

## Preprocessing

**Get the names right**

``` r
fam = paste0(sample_basename,".fam")
famtab = data.frame(fread(fam),stringsAsFactors = F)

newfam = data.frame(oldFID = famtab$V1, oldIID = famtab$V2, FID = paste0("USO",sprintf("%06d", c(1:nrow(famtab)))), IID = paste0("USO",sprintf("%06d", c(1:nrow(famtab))),1), PHENO = 1, stringsAsFactors = F)

#write.table(newfam, paste0(dir,"/USociety_IDS.txt"), quote = F, row.names = F, col.names = F)
```

awk ‘{print $1,$2,$3,$4}’ Input/USociety\_final/USociety\_IDS.txt \>
tmpIDS.txt; ./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety –update-ids tmpIDS.txt –make-bed –out
TEMP1; awk ‘{print $3,$4,$5}’ Input/USociety\_final/USociety\_IDS.txt \>
tmpPHENO.txt; ./resources/software/plink1.9/plink –bfile TEMP1 –pheno
tmpPHENO.txt –make-bed –out
Input/USociety\_final/USociety\_updated\_all; rm temp\* TEMP\*

**Here we need to specify some locations for the analysis again**

``` r
sample_basename = "Input/USociety_final/USociety_updated_all"

dir = dirname(sample_basename)

qcdir = file.path(dir, "Qc/")

finalname = "USociety_clean"
```

## Quality Control

Run PLINK Linux commands to:

**1)Generate the missing rate/calling rate, heterozygocity and sex
checking files **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_updated\_all –missing –het –check-sex
–out Input/USociety\_final/Qc/USociety\_updated\_all

**2)Prune SNPs in relative LD so we can correctly calculate IBD **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_updated\_all –indep-pairwise 50 5 0.5
–out Input/USociety\_final/Qc/USociety\_updated\_all

**3)IBD on pruned variants **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_updated\_all –maf 0.05 –genome –out
Input/USociety\_final/Qc/USociety\_updated\_all –extract
Input/USociety\_final/Qc/USociety\_updated\_all.prune.in

### Input and Read data

Import all the files that plink generated into R. Then we can plot and
do some statistical analyses on them.

``` r
#Define the file names
imiss = file.path(qcdir,paste0(basename(sample_basename),'.imiss'))
het = file.path(qcdir,paste0(basename(sample_basename),'.het'))
sexcheck = file.path(qcdir,paste0(basename(sample_basename),'.sexcheck'))
ibd = file.path(qcdir,paste0(basename(sample_basename),'.genome'))

# Stop the script if the data file cannot be found in the current working directory:
stopifnot(file.exists(imiss))
stopifnot(file.exists(het))
stopifnot(file.exists(sexcheck))
stopifnot(file.exists(ibd))

#Load the actual data into R
miss_tab <- read.table(file=imiss, header=T, stringsAsFactors = FALSE)
het_tab <- read.table(file=het, header=T, stringsAsFactors = FALSE)
sexchk_tab <- read.table(file=sexcheck, header=T, stringsAsFactors = FALSE)
ibd_tab <- data.frame(fread(ibd, header=T), stringsAsFactors = FALSE)
```

### Calculating Heterozygocity and Calling Rate per sample

Autosomal Heterozygocity = (Number of autosomal genotypes - Observed
number of homozygotes)/Number of autosomal genotypes  
Calling Rate = 1-(Number of missing genotype calls/ Number of valid
calls)  
Outliers: Samples with Heterozygocity \< mean-3SD heterozygocity calling
rate  
Samples with Heterozygocity \> mean+3SD heterozygocity calling rate  
Samples with Calling rate \< 0.97

``` r
#Calculate the heterozygocity rate
het_rate = (het_tab$N.NM. - het_tab$O.HOM.)/het_tab$N.NM.

#Calculate the genotype calling rate
cal_rate = 1-miss_tab$F_MISS

names(het_rate) = het_tab$IID
names(cal_rate) = miss_tab$IID

#Name the outliers
outliers = unique(c(names(cal_rate[cal_rate < 0.97]), names(het_rate[het_rate > mean(het_rate)+(3*sd(het_rate)) |                het_rate <  mean(het_rate)-(3*sd(het_rate))])))
outliers_rates = data.frame(Name=outliers, Call=cal_rate[outliers], Het=het_rate[outliers])

bline_outliers = unique(c(names(cal_rate[cal_rate < 0.975]), names(het_rate[het_rate > mean(het_rate)+(2*sd(het_rate))                  | het_rate < mean(het_rate)-(2*sd(het_rate))])))
bline_outliers_rates = data.frame(Name=bline_outliers, Call=cal_rate[bline_outliers], Het=het_rate[bline_outliers])

#Draw the heterozygocity to calling rate plot
hc_plot = ggplot(miss_tab, aes(x=het_rate, y=cal_rate, label=IID)) + 
          geom_point(size=2, shape=21) +
          labs(x="Autosomal Heterozygocity Rate", y = "Calling Rate") + 
          theme_classic() +
          geom_vline(xintercept=0.2275054, colour='green', linetype = "longdash", size=0.2) +
          geom_vline(xintercept=0.1960522, colour='green', linetype = "longdash", size=0.2) +
          geom_hline(yintercept=0.97, colour='gray', linetype = "longdash", size=0.2)


outliers = outliers
hc_plot
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-157-1.png)<!-- -->
Grey horizontal line: Calling Rate=0.97 threshold Green vertical line:
±3 SD deviation from the heterozygosity rate mean of the samples
threshold 77 Outliers (will be excluded)

``` r
ibd_tab = ibd_tab[!(ibd_tab$IID1 %in% outliers | ibd_tab$IID2 %in% outliers),]
sexchk_tab = sexchk_tab[!sexchk_tab$IID %in% outliers,]
snpsex = sexchk_tab$SNPSEX
names(snpsex) = sexchk_tab$IID
```

### Check relatedness from the genotype calls by measuring Identity-by-descent (IBD)

Here we check the IBD stats of each sample pair in our samples
(Excluding the previously pointed outliers).  
In our data-set there are some related samples which we need to exclude.
We therefore expect a PI\_HAT \< 0.05 for most of the samples but not
for these related ones.

To check the relatedness between the samples we’re going to inspect the
PI\_HAT distribution:

``` r
#Plot the IBD of all samples without the previously excluded ones
#barplot(ibd_tab$PI_HAT, main="IBD plot before relateds exclusion", xlab = 'Samples', ylab='Identity-by-descent (IBD)', ylim=c(0,1))
#abline(h = 0.05, lty = 2, col='red')
#write.table(ibdth,"tempibd.txt", quote = F, row.names = F)

#IBD stats before excluding relateds:
summary(ibd_tab$PI_HAT)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.000000 0.000000 0.000000 0.003524 0.006300 1.000000

For unrelated samples the highest acceptable IBD is 0.05 (red line). The
related samples have an IBD of around 0.5 as it was expected (we knew
they were related, so this is a verification of our analysis). These
samples will be excluded. Related samples Exclusion criteria: 1) Sex:
Female \> Male 2) Genotyping Calling Rate

After relateds exclusion:

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.000000 0.000000 0.000000 0.003524 0.006300 1.000000

    ## 27  related controls will be excluded

``` r
#In the previously defined outlier list, add the related samples
#outliers = unique(c(outliers, unique(ibd_exc)))
```

Perfect\!

Create a file of these IIDs and their FIDs which will be excluded. This
file will be used with the –remove plink command.

``` r
plink_exclude = miss_tab[miss_tab$IID %in% outliers,c('FID', 'IID')]
if(length(outliers)!=nrow(plink_exclude)){stop()}
#write.table(plink_exclude, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE)

curfam = read.table(paste0(sample_basename,".fam"),stringsAsFactors = F)
curfam$newsex = snpsex[curfam$V2]
#write.table(curfam[,c(1,2,3,4,7,6)], paste0(sample_basename,".fam"), sep=' ', row.names = FALSE, quote=FALSE, col.names = FALSE)
```

## PCA Analysis

### Create new plink clean set of files without the excluded samples:

We need to remove all samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_updated\_all –remove tempexcIDS.txt
–filter-females –make-bed –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd; rm temp*;
rm toexc*;

### LD prune and remove any missing calls from the new clean data-set (Like we did in the beginning of the analysis)

***1)Prune SNPs in relative LD and remove missing calls so we can
correctly calculate Principal Components later ***

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd
–indep-pairwise 50 5 0.5 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldprune

**2) Remove the pruned variants from the dataset **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd –maf 0.01
–geno 0.05 –genome –extract
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldprune.prune.in
–make-bed -out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldpruned

### Calculation of principal components (PCA Plot) of our clean and pruned Data-Set

We will generate the PCA plot to check how the samples are clustered
based on their genotype data. Any batch effects such as ethnicity will
be revealed. We’ll be using GCTA software for this (check Dependencies.R
script).

**First, we need to convert our plink files to GCTA grm format**

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldpruned
–make-grm –thread-num 10 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldpruned

**Now we can use the grm files we created to perform the GCTA PC
analysis. The output eigenvec and eigenval files are the ones we need
for the PCA plot**

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldpruned
–pca 4 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd.ldpruned

\#\#\#PCA Plot on the LD pruned data-set

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)
pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)
pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))
colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
pcs$group = NA
pcs[startsWith(pcs$FID,"U"),]$group = "USociety"


#PC1 vs PC2
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC2") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          theme_classic()

#PC2 vs PC3
pca2 = ggplot(pcs, aes(x=PC2, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC2", y = "PC3") + #geom_text(aes(label=ifelse((abs(PC2)>4*IQR(PC2)|abs(PC3)>4*IQR(PC3)),IID,"")), vjust=-0.6) +
          theme_classic()

#PC1 vs PC3
pca3 = ggplot(pcs, aes(x=PC1, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC3") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC3)>4*IQR(PC3)),IID,"")), vjust=-0.6) +
          theme_classic()

#Table with EigenValues and all the PCA plots
pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues') %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Eigenvalues

</caption>

<thead>

<tr>

<th style="text-align:right;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

5.18538

</td>

</tr>

<tr>

<td style="text-align:right;">

3.47582

</td>

</tr>

<tr>

<td style="text-align:right;">

3.04946

</td>

</tr>

<tr>

<td style="text-align:right;">

2.67612

</td>

</tr>

<tr>

<td style="text-align:right;">

2.58114

</td>

</tr>

</tbody>

</table>

``` r
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-180-1.png)<!-- -->

``` r
pca2
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-180-2.png)<!-- -->

``` r
pca3
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-180-3.png)<!-- -->

### Merge our data with three Hapmap populations; CEU; Yoruban and CHB.

**PCA analysis of GWAS & HapMap samples** We want to only select the
samples which cluster together with the CEU hapmap population.
**1)Extract HapMap SNPs **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd –extract
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.no-at-cg-snps.txt
–make-bed –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps

**2)We will merge the HapMap and Batch2 datasets but multi-allelic SNPs
will be created: If an allele is the major in one dataset and the minor
in the other dataset, then plink marks them as multi-allelic SNPs which
we will remove.**

**Initially merge the HapMap and Batch2 datasets:**

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed

**Remove the multi-allelic SNPs.**

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps
–exclude plink.missnp –make-bed –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed
–out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Prune SNPs in relative LD and remove missing calls so we can correctly
calculate Principal Components later**

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA
–indep-pairwise 50 5 0.5 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Remove the pruned variants from the dataset **

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA
–maf 0.01 –geno 0.05 –extract
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.prune.in
–make-bed -out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Convert our plink files to GCTA grm format **

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–make-grm –autosome –thread-num 10 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Now we can use the grm files we created to perform the GCTA PC
analysis. The output eigenvec and eigenval files are the ones we need
for the PCA plot **

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–pca 4 –out
Input/USociety\_final/Qc/USociety\_updated\_all.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)
pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)

colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))


pcs$group = NA
pcs[startsWith(pcs$FID,"U"),]$group = "USociety"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"1"),]$group = "CEU"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"NA"),]$group = "CHB"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"Y"),]$group = 'YRI'
pcs$group = factor(pcs$group)


pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))


#pc_vals = read.table(paste0(plink_dir,out_suffix,'_hapmap-snps2_for_PCA_pruned.eigenval'), header=FALSE, stringsAsFactors = FALSE)

#colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
#colnames(pc_vals)=c('EigenValues')

#PC1 vs PC2
#pcs = pcs[pcs$group %in% c("Lipoedema_Batch1","Lipoedema_Batch2"),]
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group, label=paste0(FID,"-",IID))) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC2") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          theme_classic()
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-205-1.png)<!-- -->

``` r
pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues')  %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Eigenvalues

</caption>

<thead>

<tr>

<th style="text-align:right;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

91.45710

</td>

</tr>

<tr>

<td style="text-align:right;">

38.67590

</td>

</tr>

<tr>

<td style="text-align:right;">

3.37738

</td>

</tr>

<tr>

<td style="text-align:right;">

3.27517

</td>

</tr>

<tr>

<td style="text-align:right;">

2.83262

</td>

</tr>

</tbody>

</table>

``` r
#gg = ggplotly(pca1)
#gg
```

These are outliers around the European population and therefeore we
excluded them.

PCA After exclusion:

``` r
pcs = pcs[!pcs$IID %in% pca_outs[,2],]

#PC1 vs PC2
#pcs = pcs[pcs$group %in% c("Lipoedema_Batch1","Lipoedema_Batch2"),]
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group, label=paste0(FID,"-",IID))) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC2") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          theme_classic()
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-207-1.png)<!-- -->

``` r
pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues')  %>% kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Eigenvalues

</caption>

<thead>

<tr>

<th style="text-align:right;">

x

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

91.45710

</td>

</tr>

<tr>

<td style="text-align:right;">

38.67590

</td>

</tr>

<tr>

<td style="text-align:right;">

3.37738

</td>

</tr>

<tr>

<td style="text-align:right;">

3.27517

</td>

</tr>

<tr>

<td style="text-align:right;">

2.83262

</td>

</tr>

</tbody>

</table>

``` r
#gg = ggplotly(pca1)
```

### Create new plink clean set of files without the excluded samples and remove males:

Add all under-exclusion samples together:

``` r
final_exclude = rbind(plink_exclude,pca_outs)

#write.table(final_exclude, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE)
```

We need to remove all four samples we decided to exclude from the
analysis:

./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_updated\_all –remove tempexcIDS.txt –maf
0.01 –geno 0.05 –filter-females –hwe 1e-6 midp –make-bed –out
Input/USociety\_final/USociety\_clean; rm temp\*

# Lipoedema vs Understanding Society GWAS

Merge lipoedema cases and controls genotyping data and perform the
association analysis.

## Files management

We need to specify the locations of the ouput directory (Batch
directory) and the locations of the two plink filesets we are going to
merge

``` r
##Where are we going to run our analysis?
b1 = "Merged/Lipoedema_Batches_final/lipoedema_clean"
b2 = "Input/USociety_final/USociety_clean"

sample_basename = "Merged/USociety_vs_Lipoedema_final/lipoUSociety_merged.no_at-cg.no_inc"

dir = dirname(sample_basename)
qcdir = file.path(dir, "Qc/")
analysis_dir = file.path(dir, "Association_analysis")
meta_dir = file.path(dir, "Meta-analysis")

finalname = "lipoUSo_clean"

masterfile = "./Input/Sample_info/Glens recruits-merged_v3.csv"
lmaster = read.csv(masterfile, stringsAsFactors = F)
diagnoses <- lmaster$DIAGNOSIS.CONF.
names(diagnoses) <- lmaster$D.Number

dir.create(dir, showWarnings = FALSE)
dir.create(qcdir, showWarnings = FALSE)
dir.create(analysis_dir, showWarnings = FALSE)
dir.create(meta_dir)
```

    ## Warning in dir.create(meta_dir): 'Merged\USociety_vs_Lipoedema_final\Meta-
    ## analysis' already exists

## Load the bim and fam files of the plink filesets we want to merge

``` r
b1_bim = data.frame(fread(paste0(b1,".bim")), stringsAsFactors = FALSE)
b2_bim = data.frame(fread(paste0(b2,".bim")), stringsAsFactors = FALSE)

b1_fam = data.frame(fread(paste0(b1,".fam")), stringsAsFactors = FALSE)
b2_fam = data.frame(fread(paste0(b2,".fam")), stringsAsFactors = FALSE)

rownames(b1_bim)=b1_bim$V2
rownames(b2_bim)=b2_bim$V2

#Add a unique genomic position identifier (pid) with format [no of chromosome]chr[genomic coordinate]
b1_bim$pid = paste0(b1_bim$V1,'chr',b1_bim$V4)
b2_bim$pid = paste0(b2_bim$V1,'chr',b2_bim$V4)
```

## Preprocessing

We need to only keep the common variants between the two batches.
Variants not covered in both filesets will generate unwanted biases and
batch effects.

**Keep in Batch1 only the variants which are shared between Batch2 and
Batch1**

less Input/USociety\_final/USociety\_clean.bim | awk ‘{print $2}’ \>
tempIDS.txt; ./resources/software/plink1.9/plink –bfile
Merged/Lipoedema\_Batches\_final/lipoedema\_clean –extract tempIDS.txt
–make-bed –out Merged/USociety\_vs\_Lipoedema\_final/b1.shared

**Keep in Batch2 only the variants which are shared between Batch2 and
Batch1**

less Merged/Lipoedema\_Batches\_final/lipoedema\_clean.bim | awk ‘{print
$2}’ \> tempIDS.txt; ./resources/software/plink1.9/plink –bfile
Input/USociety\_final/USociety\_clean –filter-females –extract
tempIDS.txt –make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/b2.shared

## Merging the Lipoedema and Control Understanding Society filesets together

**1)Merge Batch2 and Filtered Batch1 (This step will fail due to strand
flipping issues but it will create the usefull missnp file which will be
used for strand flipping.)**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/b2.shared –bmerge
Merged/USociety\_vs\_Lipoedema\_final/b1.shared

**2)Use the missnp file to filter flip the problematic snps in batch1
and then redo the merging.**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/b1.shared –flip plink.missnp
–make-bed –out Merged/USociety\_vs\_Lipoedema\_final/b1.shared.flipped

Now we can actually merge batch2 with the flipped version of batch1
**3)Merge the flipped batch1 with the batch2 data-set **  
In this step **WE ALSO REMOVE NON ACGT SNPS from the merged dataset**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/b2.shared –bmerge
Merged/USociety\_vs\_Lipoedema\_final/b1.shared.flipped –merge-mode 1
–snps-only just-acgt –make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/b1b2merged

**Merge is now completed.**

## Remove T -\> A, A -\> T, G -\> C, C -\> G SNPs

T -\> A, A -\> T, G -\> C, C -\> G SNPs are particularly hard to resolve
when merging datasets. In this case we are going to remove them from the
merged dataset.

``` r
cat(paste0(length(unique(at_cg))," T -> A, A -> T, G -> C, C -> G SNPs are being removed:\n",paste(unique(c(at_cg)),collapse=", ")))
```

    ## 0 T -> A, A -> T, G -> C, C -> G SNPs are being removed:

Use plink to remove these SNPs

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/b1b2merged –allow-no-sex –exclude
tempIDS.txt –make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/lipoUSociety\_merged.no\_at-cg.no\_inc

## Quality Control

**1)Generate the missing rate/calling rate, heterozygocity and sex
checking files **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/lipoUSociety\_merged.no\_at-cg.no\_inc
–missing –het –check-sex –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc

**2)Prune SNPs in relative LD so we can correctly calculate IBD **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/lipoUSociety\_merged.no\_at-cg.no\_inc
–indep-pairwise 50 5 0.5 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc

**3)IBD on pruned variants **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/lipoUSociety\_merged.no\_at-cg.no\_inc
–maf 0.05 –genome –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc
–extract
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.prune.in

## Input and Read data

Import all the files that plink generated into R. Then we can plot and
do some statistical analyses on them.

``` r
#Define the file names
imiss = file.path(qcdir,paste0(basename(sample_basename),'.imiss'))
het = file.path(qcdir,paste0(basename(sample_basename),'.het'))
sexcheck = file.path(qcdir,paste0(basename(sample_basename),'.sexcheck'))
ibd = file.path(qcdir,paste0(basename(sample_basename),'.genome'))

# Stop the script if the data file cannot be found in the current working directory:
stopifnot(file.exists(imiss))
stopifnot(file.exists(het))
stopifnot(file.exists(sexcheck))
stopifnot(file.exists(ibd))

#Load the actual data into R
miss_tab <- read.table(file=imiss, header=T, stringsAsFactors = FALSE)
het_tab <- read.table(file=het, header=T, stringsAsFactors = FALSE)
sexchk_tab <- read.table(file=sexcheck, header=T, stringsAsFactors = FALSE)
ibd_tab <- data.frame(fread(file=ibd, header=T), stringsAsFactors = FALSE)
```

# QC Measurements

## Calculating Heterozygocity and Calling Rate per sample

Autosomal Heterozygocity = (Number of autosomal genotypes - Observed
number of homozygotes)/Number of autosomal genotypes  
Calling Rate = 1-(Number of missing genotype calls/ Number of valid
calls)  
Outliers: Samples with Heterozygocity \< mean-3SD heterozygocity calling
rate  
Samples with Heterozygocity \> mean+3SD heterozygocity calling rate  
Samples with Calling rate \< 0.97

``` r
#Calculate the heterozygocity rate
het_rate = (het_tab$N.NM. - het_tab$O.HOM.)/het_tab$N.NM.

#Calculate the genotype calling rate
cal_rate = 1-miss_tab$F_MISS

names(het_rate) = het_tab$IID
names(cal_rate) = miss_tab$IID

#Name the outliers
outliers = unique(c(names(cal_rate[cal_rate < 0.97]), names(het_rate[het_rate > mean(het_rate)+(3*sd(het_rate)) |                het_rate <  mean(het_rate)-(3*sd(het_rate))])))
outliers_rates = data.frame(Name=outliers, Call=cal_rate[outliers], Het=het_rate[outliers])

bline_outliers = unique(c(names(cal_rate[cal_rate < 0.975]), names(het_rate[het_rate > mean(het_rate)+(2*sd(het_rate))                  | het_rate < mean(het_rate)-(2*sd(het_rate))])))
bline_outliers_rates = data.frame(Name=bline_outliers, Call=cal_rate[bline_outliers], Het=het_rate[bline_outliers])

#Draw the heterozygocity to calling rate plot
hc_plot = ggplot(miss_tab, aes(x=het_rate, y=cal_rate)) + 
          geom_point(size=2, shape=21) +
          labs(x="Autosomal Heterozygocity Rate", y = "Calling Rate") + 
          theme_classic() +
          geom_vline(xintercept=mean(het_rate)-(2*sd(het_rate)), colour='blue', linetype = "longdash", size=0.2) +
          geom_vline(xintercept=mean(het_rate)+(2*sd(het_rate)), colour='blue', linetype = "longdash", size=0.2)+
          geom_hline(yintercept=0.97, colour='gray', linetype = "longdash", size=0.2) #+
          #geom_text(data=outliers_rates, aes(x=Het, y=Call, label=Name), hjust=0, vjust=+1.5, size=3)

outliers = outliers

hc_plot
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-244-1.png)<!-- -->

    ## Grey horizontal line: Calling Rate=0.97 threshold
    ## Green vertical line: ±3 SD deviation from the heterozygosity rate mean of the samples threshold
    ## 34 outliers will be excluded:
    ## 34 are Controls and 0 are Cases.

``` r
plink_exclude1 = miss_tab[miss_tab$IID %in% outliers,c('FID', 'IID')]

het_tab = het_tab[!(het_tab$IID %in% outliers),]
miss_tab = miss_tab[!(miss_tab$IID %in% outliers),]
ibd_tab = ibd_tab[!(ibd_tab$IID1 %in% outliers | ibd_tab$IID2 %in% outliers),]
sexchk_tab = sexchk_tab[!sexchk_tab$IID %in% outliers,]
```

Create a file of these IIDs and their FIDs which will be excluded. This
file will be used with the –remove plink command.

``` r
if(length(outliers)!=nrow(plink_exclude1)){stop()}
#write.table(plink_exclude1, 'tempexcIDS.txt', sep='\t', row.names = FALSE, quote=FALSE, col.names = F)
#write.table(plink_exclude1, file.path(qcdir,'toremovesamples.txt'), sep='\t', row.names = FALSE, quote=FALSE, col.names = F)
```

## Check Sex

Sex assignments have been checked before.

## IBD check

We have removed all the relateds Lipoedema and Control USociety samples
so we don’t expect any more relateds in the cohort.

Let’s check the IBD summary statistics (PI\_HAT \< 0.05).

``` r
summary(ibd_tab$PI_HAT)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.000000 0.000000 0.000000 0.003466 0.006100 0.048200

Maximum PI\_HAT is \< 0.05 so there are no related people in the
data-set.

## Perform PCA

### Create new plink clean set of files without the excluded samples:

We need to remove all samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/lipoUSociety\_merged.no\_at-cg.no\_inc
–remove Merged/USociety\_vs\_Lipoedema\_final/Qc/toremovesamples.txt
–make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd;
rm temp*; rm toexc*;

### LD prune and remove any missing calls from the new clean data-set (Like we did in the beginning of the analysis)

***1)Prune SNPs in relative LD and remove missing calls so we can
correctly calculate Principal Components later ***

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd
–indep-pairwise 50 5 0.5 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldprune

**2) Remove the pruned variants from the dataset **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd
–maf 0.01 –geno 0.05 –extract
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldprune.prune.in
–make-bed -out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldpruned

## Calculation of principal components (PCA Plot) of our clean and pruned Data-Set

We will generate the PCA plot to check how the samples are clustered
based on their genotype data. Any batch effects such as ethnicity will
be revealed. We’ll be using GCTA software for this (check Dependencies.R
script).

First, we need to convert our plink files to GCTA grm format

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldpruned
–make-grm –thread-num 10 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldpruned

Now we can use the grm files we created to perform the GCTA PC analysis.
The output eigenvec and eigenval files are the ones we need for the PCA
plot

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldpruned
–pca 8 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd.ldpruned

### PCA Plot on the LD pruned data-set

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)

#write.table(pcs[,c(1:4)],file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.2pcs.eigenvec")), quote=F, col.names = F, row.names = F)
#write.table(pcs[,c(1:5)],file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.3pcs.eigenvec")), quote=F, col.names = F, row.names = F)
#write.table(pcs[,c(1:6)],file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.4pcs.eigenvec")), quote=F, col.names = F, row.names = F)

pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd.ldpruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)
pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))
colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
pcs$group = NA
pcs[startsWith(pcs$FID,"B1_"),]$group = "Lipoedema_B1"
pcs[startsWith(pcs$FID,"B2_"),]$group = "Lipoedema_B2"
pcs[startsWith(pcs$FID,'USO'),]$group = "USociety"

#PC1 vs PC2
pca1 = ggplot(pcs, aes(x=PC1, y=PC2, color=group)) + 
          geom_point(size=2, shape=21) + geom_vline(xintercept = -0.1) +
          labs(x="PC1", y = "PC2") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          theme_classic()

#PC2 vs PC3
pca2 = ggplot(pcs, aes(x=PC2, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC2", y = "PC3") + #geom_text(aes(label=ifelse((abs(PC2)>4*IQR(PC2)|abs(PC3)>4*IQR(PC3)),IID,"")), vjust=-0.6) +
          theme_classic()

#PC1 vs PC3
pca3 = ggplot(pcs, aes(x=PC1, y=PC3, color=group)) + 
          geom_point(size=2, shape=21) +
          labs(x="PC1", y = "PC3") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC3)>4*IQR(PC3)),IID,"")), vjust=-0.6) +
          theme_classic()

#Table with EigenValues and all the PCA plots
#pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues') %>% kable_styling()
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-264-1.png)<!-- -->

``` r
pca2
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-264-2.png)<!-- -->

``` r
pca3
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-264-3.png)<!-- -->

``` r
pcouts1 = pcs[pcs$PC1 < (-0.11),]$IID
```

Since we have excluded all PCA outliers in batch-specific analyeses, as
expected we have no outliers in our PCA plots which need to be excluded.
However we are going to run the PCA analysis again using the HapMap
samples as population references to see how the analysis’ samples
cluster with the hapmap populations.

## PCA analysis of GWAS & HapMap samples

We expect all our samples to cluster together with the Central European
samples.

**1)Extract HapMap SNPs **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd
–extract ./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.no-at-cg-snps.txt
–make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps

**2)We will merge the HapMap and Batch2 datasets but multi-allelic SNPs
will be created: If an allele is the major in one dataset and the minor
in the other dataset, then plink marks them as multi-allelic SNPs which
we will remove. ** **-Initially merge the HapMap and Batch2 datasets: **

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed

**Remove the multi-allelic SNPs**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps
–exclude plink.missnp –make-bed –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2
–bmerge
./resources/hapmap/hapmap3r2\_CEU.CHB.JPT.YRI.founders.no-at-cg-snps\_hg19\_bad\_snps\_removed
–out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Prune SNPs in relative LD and remove missing calls so we can correctly
calculate Principal Components later**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA
–indep-pairwise 50 5 0.5 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA

**Remove the pruned variants from the dataset**

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA
–maf 0.01 –geno 0.05 –extract
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.prune.in
–make-bed -out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Convert our plink files to GCTA grm format **

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–make-grm –autosome –thread-num 10 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

**Now we can use the grm files we created to perform the GCTA PC
analysis. The output eigenvec and eigenval files are the ones we need
for the PCA plot **

./resources/software/gcta\_1.93.2beta/gcta64 –grm
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned
–pca 4 –out
Merged/USociety\_vs\_Lipoedema\_final/Qc/lipoUSociety\_merged.no\_at-cg.no\_inc.cal.het.ibd\_hapmap-snps2\_for\_PCA.pruned

``` r
#Load the eigenvectors and eigenvalues to draw the plot
pcs = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenvec")), header=FALSE, stringsAsFactors = FALSE)
pc_vals = read.table(file.path(qcdir,paste0(basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.eigenval")), header=FALSE, stringsAsFactors = FALSE)
colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
pcs$group = NA
pcs[startsWith(pcs$FID,"B1_"),]$group = "Lipoedema"
pcs[startsWith(pcs$FID,"B2_"),]$group = "Lipoedema"
pcs[startsWith(pcs$FID,"USO"),]$group = "Control"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"1"),]$group = "CEU"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"NA"),]$group = "CHB"
pcs[startsWith(pcs$IID,"NA") & startsWith(pcs$FID,"Y"),]$group = 'YRI'
#
pcs2 = rbind(pcs[pcs$group=="Control",],pcs[pcs$group=="Lipoedema",],pcs[pcs$group=="CEU",],pcs[pcs$group=="CHB",],pcs[pcs$group=="YRI",])

pc_vals_tab = data.frame(PC_Eigenvalues=pc_vals[c(1:5),])
rownames(pc_vals_tab)=paste0('PC', c(1:nrow(pc_vals_tab)))

pcs2$alpha = 0.2
pcs2[pcs2$group=="Lipoedema",]$alpha = 0.5
pcs2[pcs2$group=="Control",]$alpha = 0.3
pcs2[pcs2$group=="YRI",]$alpha = 0.3
pcs2[pcs2$group=="CHB",]$alpha = 0.3
pcs2$group = factor(pcs2$group, levels = c("Control","Lipoedema","CEU","CHB","YRI"))
 

#pc_vals = read.table(paste0(plink_dir,out_suffix,'_hapmap-snps2_for_PCA_pruned.eigenval'), header=FALSE, stringsAsFactors = FALSE)

#colnames(pcs)=c('FID','IID',paste0('PC', c(1:(ncol(pcs)-2))))
#colnames(pc_vals)=c('EigenValues')

#PC1 vs PC2
#pcs = pcs[pcs$group %in% c("Lipoedema_Batch1","Lipoedema_Batch2"),]
pca1 = ggplot(pcs2, aes(x=PC1, y=PC2, color=group, label=IID)) + 
          geom_point(size=5, shape=19, aes(alpha=alpha)) +
          labs(x=paste0("PC1: ",round(pc_vals$V1[1],1),"% Variance"), y = paste0("PC2: ",round(pc_vals$V1[2],1),"% Variance")) + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
          scale_colour_manual(values=c("#56B4E9","#E69F00","#FF0606","#009E73","#F0E442")) +
          guides(fill=guide_legend(title="Cohort")) +
          theme_bw()
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-289-1.png)<!-- -->

``` r
#pc_vals_tab[c(1:5),] %>% kable(caption = 'Eigenvalues')  %>% kable_styling()
#gg = ggplotly(pca1)
#gg

#tiff(paste0(qcdir,basename(sample_basename),".cal.het.ibd_hapmap-snps2_for_PCA.pruned.pca_plot2.tiff"),width=2600,height=1500, res = 600, pointsize = 12)
pca1
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-289-2.png)<!-- -->

``` r
dev.off()
```

    ## null device 
    ##           1

``` r
plink_exclude_all = plink_exclude1
#write.table(plink_exclude_all, file.path(qcdir,'toremovesamples_all_final.txt'), sep='\t', row.names = FALSE, quote=FALSE, col.names = F)
```

\#\#Create the final clean set of files: We need to remove all the
samples we decided to exclude from the analysis:

./resources/software/plink1.9/plink –bfile
Merged/USociety\_vs\_Lipoedema\_final/b1b2merged –geno 0.05 –maf 0.01
–hwe 1e-6 midp –exclude
Merged/USociety\_vs\_Lipoedema\_final/Qc/atcgIDS.txt –remove
Merged/USociety\_vs\_Lipoedema\_final/Qc/toremovesamples\_all\_final.txt
–make-bed –out Merged/USociety\_vs\_Lipoedema\_final/lipoUSo\_clean

# SNP-Based Heritability estimation

### Files management

We need to specify the locations of the ouput directory (Batch
directory) and the locations of the two plink filesets we are going to
merge

``` r
##Where are we going to run our analysis?

sample_basename = "./Merged/USociety_vs_Lipoedema_final/lipoUSo_clean"

dir = dirname(sample_basename)
herdir = file.path(dir, "Heritability_estimation/")
qcdir = file.path(dir, paste0(herdir,"/Qc/"))

finalname = "lipoUSo_clean"

dir.create(dir, showWarnings = FALSE)
dir.create(qcdir, showWarnings = FALSE)
dir.create(herdir, showWarnings = FALSE)
```

**Load the fam file of the plink filesets we want to perform REML on**

``` r
ffam = data.frame(fread(paste0(sample_basename,".fam")), stringsAsFactors = FALSE)

gcta_pheno = cbind(ffam[,c("V1","V2")],ffam$V6-1)
write.table(gcta_pheno, paste0(herdir,"/gcta_phenos.txt"), sep="\t", row.names=F, col.names=F, quote=F)
```

### GCTA analysis

1)First, we need to convert our plink files to GCTA grm format

./resources/software/gcta\_1.93.2beta/gcta64 –bfile
./Merged/USociety\_vs\_Lipoedema\_final/lipoUSo\_clean –make-grm
–thread-num 10 –out
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/lipoUSo\_clean

2)  GCTA-GREML analysis: estimating the variance explained by the SNPs

**- Using prevalence=0.05**

./resources/software/gcta\_1.93.2beta/gcta64 –grm
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/lipoUSo\_clean
–pheno
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/gcta\_phenos.txt
–reml –prevalence 0.05 –out
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/lipoUSo\_clean.prev005

**- Using prevalence=0.1**

./resources/software/gcta\_1.93.2beta/gcta64 –grm
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/lipoUSo\_clean
–pheno
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/gcta\_phenos.txt
–reml –prevalence 0.1 –out
./Merged/USociety\_vs\_Lipoedema\_final/Heritability\_estimation/lipoUSo\_clean.prev01

# Association Analysis

### Loading necessary functions to parse the association analysis results

``` r
results_annotation <- function(res_path_glm, res_path_19) {
  path_prefix = tools::file_path_sans_ext(res_path_glm)
  results2 <- data.frame(fread(res_path_glm, header=T), stringsAsFactors = FALSE)
  results3 <- data.frame(fread(res_path_19, header=T), stringsAsFactors = FALSE)
  
  results2 <- merge(results2, results3[,c("SNP","P","OR","SE","L95","U95")], by.x="ID", by.y="SNP", all.x=T, all.y =F,
                    suffixes = c("_glm","_assoc"))
  
  if("MT" %in% results2$X.CHROM){
    results2[results2$X.CHROM == "MT",]$X.CHROM = 25
  }
  
  if("X" %in% results2$X.CHROM){
    results2[results2$X.CHROM == "X",]$X.CHROM = 23
  }
  
  if("Y" %in% results2$X.CHROM){
    results2[results2$X.CHROM == "Y",]$X.CHROM = 24
  }
  
  results2$X.CHROM = as.numeric(results2$X.CHROM)
  
  results2$nexus = "dbsnp"
  write.table(head(results2[order(results2$P_assoc),c("nexus","ID")],300), paste0(path_prefix,".for_nexus.txt"),  sep="\t", quote=F,
              col.names = F, row.names = F)
  
  readline(prompt="Upload the for nexus file on nexus. save the annotation and then Press [enter] to continue")
  
  results2_genemapping = read.table(paste0(path_prefix,".gene_mapping.txt"),stringsAsFactors = F,sep="\t",header=T)

  results2_genemapping$Gene = ifelse(results2_genemapping$Overlapped.Gene=="None",paste0(results2_genemapping$Nearest.Upstream.Gene,";",results2_genemapping$Nearest.Downstream.Gene),results2_genemapping$Overlapped.Gene)
  
  results2_top300 = merge(results2, results2_genemapping[,c("Variation.ID","Gene","Annotation")],by.x="ID",by.y="Variation.ID", all.x=F)
  
  write.table(results2_top300, paste0(path_prefix,".top300annot.txt"),  sep="\t", quote=F, col.names = F, row.names = F)
  
  return(results2_top300)
}

make_plots_ready_plink2 <- function(res_path){
  
  path_prefix = tools::file_path_sans_ext(res_path)  
  results2 <- data.frame(fread(res_path, header=T), stringsAsFactors = FALSE)
  
  if("X.CHROM" %in% colnames(results2)){
    results2 = results2[!results2$CHROM %in% c("MT","X"),]
    results2$X.CHROM = as.numeric(results2$X.CHROM)
  }
  if("CHR" %in% colnames(results2)){
    results2 = results2[!results2$CHR %in% c("MT","X"),]
    results2$CHR = as.numeric(results2$CHR)
  }  
  
  fresults=results2
  
  return(fresults)
}
```

### Running the Association Analysis between Cases and Controls

We are going to run both the simple logistic regression analysis on
PLINK 1.9 which will provide the OR and the PValue we are counting on
but we are also goint to run the PLINK 2.0 glm analysis without
covariates to obtain some nice stats that only PLINK 2.0 outputs:

**Run the Association Analysis**

``` r
temp_comm = paste0(plink,' --bfile ',file.path(dir,finalname),' --logistic --ci 0.95 ',' --out ',
                   file.path(analysis_dir,finalname),"; \n",
                   plink2,' --bfile ',file.path(dir,finalname),
                   ' --glm allow-no-covars cols=chrom,pos,ref,alt1,alt,a1countcc,totallelecc,gcountcc,nobs,beta,orbeta,tz,p',
                   ' --out ',file.path(analysis_dir,finalname))
```

./resources/software/plink1.9/plink –bfile
./Merged/USociety\_vs\_Lipoedema\_final/lipoUSo\_clean –logistic –ci
0.95 –out
Merged/USociety\_vs\_Lipoedema\_final/Association\_analysis/lipoUSo\_clean;
./resources/software/plink2/plink2 –bfile
./Merged/USociety\_vs\_Lipoedema\_final/lipoUSo\_clean –glm
allow-no-covars
cols=chrom,pos,ref,alt1,alt,a1countcc,totallelecc,gcountcc,nobs,beta,orbeta,tz,p
–out
Merged/USociety\_vs\_Lipoedema\_final/Association\_analysis/lipoUSo\_clean

**Calculate genomic inflation**

We’re gonna run the same plink –logistic association but also using
–adjust so the lambda inflation factor will be calculated

``` r
temp_comm = paste0(plink,' --bfile ',file.path(dir,finalname),' --adjust --logistic --ci 0.95 ',' --out ',
                   file.path(analysis_dir,finalname),'_for_lambda')
```

./resources/software/plink1.9/plink –bfile
./Merged/USociety\_vs\_Lipoedema\_final/lipoUSo\_clean –adjust –logistic
–ci 0.95 –out
Merged/USociety\_vs\_Lipoedema\_final/Association\_analysis/lipoUSo\_clean\_for\_lambda

**Association Analysis Results**

``` r
#Load the assocation results from both --logistic and --glm analyses and merge them into R
if(run==T){
  strict_top300 <- results_annotation(paste0(file.path(analysis_dir,finalname),".PHENO1.glm.logistic.hybrid"),
                                      paste0(file.path(analysis_dir,finalname),".assoc.logistic"))
  
  strict_top30_assoc <- head(strict_top300[order(strict_top300$P_assoc),],35)
  strict_top30_glm <- head(strict_top300[order(strict_top300$P_glm),],35)
  
  strict_top30 <- strict_top300[strict_top300$ID %in% unique(strict_top30_assoc$ID,strict_top30_glm$ID),]
  strict_top30 <- strict_top30[!duplicated(strict_top30$ID),]
  strict_top30 <- strict_top30[order(strict_top30$P_assoc),]
  
  #Write the top30 snps
  write.table(strict_top30, paste0(file.path(analysis_dir,finalname),".PHENO1.glm.logistic",".top30annot.txt"), 
              sep="\t", quote=F, col.names = T, row.names = F)

  #write_the_file_for the meta-analysis
  discovery_meta = cbind(strict_top30[,c(1,7,4)],strict_top30$CASE_ALLELE_CT/2,strict_top30[,c(23,24,22)])
  colnames(discovery_meta) = c("SNP","A1","REF","N","OR","SE","P")
  write.table(discovery_meta, 
              paste0(meta_dir,"/lipoUSo_clean.PHENO1.glm.logistic.top30annot.case_weight.or.metal.txt"),sep="\t",
              quote=F,row.names = F)
}

#Get data ready for plot
strict_results <- data.frame(fread(paste0(file.path(analysis_dir,finalname),".assoc.logistic"), header=T), stringsAsFactors = FALSE)
strict_results_top <- head(strict_results[order(strict_results$P),],10000)
plot_strict = strict_results[,c(2,1,3,12)]
plot_strict = plot_strict[!plot_strict$CHR=='26',]
plot_strict[plot_strict$CHR=='23',]$CHR = "X"

plot_strict_top = head(plot_strict[order(plot_strict$P),],10000)

#highlight snps with meta_p<1e-4 and same direction
highlight_snps=c("rs1409440","rs7994616","rs11616618","rs9308098","rs10521271","rs11511253","rs4554078","rs2019845","rs10499948","rs16825349")

#Draw manhattan plot
par(mar=c(5.1, 5.1, 4.1, 2.1))
CMplot(plot_strict,plot.type="m",threshold=c(1e-8,1e-4),threshold.lty=c(0,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","black"),amplify=F,
       signal.col=NULL,highlight.col=NULL,file.output = F,dpi=900,
       pch=21,cex=.9,ylim=c(0,6),col=c("#648fff","#ffb000"),signal.pch = 19,signal.cex=.9,
       highlight = highlight_snps, highlight.cex = 1.5, highlight.pch = 19, cex.axis = 1, cex.lab=1.3)
```

    ##  (warning: all phenotypes will use the same ylim.)
    ##  Rectangular-Manhattan Plotting P.

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-312-1.png)<!-- -->

``` r
#Draw qq plot
par(mar=c(5.1, 5.1, 4.1, 2.1))
CMplot(plot_strict,col = "#648fff",plot.type="q",box=FALSE,file.output = F,file='tiff',memo="",dpi=900,
conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
width=8,height=5)
```

    ##  QQ Plotting P.

``` r
text(5.3,0,"lambda = 1.00364",cex = 1.5)
```

![](lipoedema_supplementary_files/figure-gfm/unnamed-chunk-312-2.png)<!-- -->

### Meta-analysis using GEL Data

We’re gonna perform meta analysis on METAL software. Association Results
for the 25 variants crossing the suggestive cutoff from the discvoery
study and the results from the replication study for the same variants
have been reformed in a format accepted by the METAL software and the
new files have been placed in the Meta-analysis directory.

**Create the METAL script**

``` r
if(run==T){
  sink(file.path(meta_dir,"metal_analysis.or.txt"))
  cat("SCHEME STDERR")
  cat("\n\n")
  cat("MARKER\tSNP
  WEIGHT\tN
  ALLELE\tA1 REF
  EFFECT\tlog(OR)
  STDERR\tSE
  PVAL\tP\n")
  cat("\n")
  cat(paste0("PROCESS ",file.path(meta_dir,finalname),
             ".PHENO1.glm.logistic.top30annot.case_weight.or.metal.txt\n"))
  cat(paste0("PROCESS ",meta_dir,"/gel_lipo_gwas_strict_all_assoc_for_meta.txt\n"))
  cat("\n")
  cat(paste0("OUTFILE ",meta_dir,"/final_meta-analysis_ .tbl\n\n"))
  cat("ANALYZE")
  sink()
}
```

**Run METAL software**

``` r
temp_comm = paste0(metal," ",file.path(meta_dir,"metal_analysis.or.txt"))
```

./resources/software/generic-metal/metal
Merged/USociety\_vs\_Lipoedema\_final/Meta-analysis/metal\_analysis.or.txt

**Load and merge all data from discovery, replication and meta-analysis
analyses**

``` r
if(run==T){
  discovery = read.table(paste0(file.path(analysis_dir,finalname),".PHENO1.glm.logistic.top30annot.txt"),
                         sep="\t", stringsAsFactors = F, header=1)
  
  
  replication = read.table(paste0(meta_dir,"/gel_lipo_gwas_strict_all_assoc.txt"),
                           stringsAsFactors = F,sep="\t",header=1)
  meta_analysis = read.table(paste0(meta_dir,"/final_meta-analysis_1.tbl"),
                             stringsAsFactors = F, sep="\t",header=1)
  
  disc_rep = merge(discovery, replication, by.x="ID",by.y="ID",suffixes=c("_SGUL","_GEL"))
  final_results = merge(disc_rep, meta_analysis, by.x="ID", by.y="MarkerName")
  
  write.table(final_results,paste0(meta_dir,"/final_all_merged.tsv"),sep="\t",quote=F,row.names = F)
}
```

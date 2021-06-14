import os
import sys
import time
import glob
import subprocess
import re
import random
import string
import pandas as pd
import csv

indir = "/re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/agg_v2/data_requests"
outdir = indir+"/Output"
rundir = indir+"/Runs"
inpdir = indir+"/Input"

files_map = inpdir + "/" + "gel_genomes_map.txt"

project = "lipo_gwas_strict_hits"

from dependencies import *
from utils import *

threads=1

map_df = pd.read_csv(files_map, sep='\t', header=0)
#map_df.columns = ["chr","start","end","file"]
map_df = map_df.sort_values(by=['chr', 'start'])
map_df = map_df.reset_index(drop=True)



project_out,project_run,project_js,project_tmp,genomic_data,interval_data = tidydirs(project, outdir, rundir)

#genomic_vcfs = glob.glob(aggvcfs+"/*vcf.gz")
genomic_vcfs = glob.glob("/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data/*vcf.gz")
#print(genomic_vcfs)
gvcfs_dict = {x:[x.split(".")[0].split("_")[-3],x.split(".")[0].split("_")[-2],x.split(".")[0].split("_")[-1]] for x in genomic_vcfs}

snps_file = inpdir + "/lipo_gwas_strict_replication" + "/" + "gwas_top_hits.txt"
fam_file_y4 = inpdir + "/lipo_gwas_strict_replication" + "/" + "masterfile.y1y2y3y4.fam"

exclude_ancestry = inpdir + "/lipo_gwas_strict_replication" + "/" + "ancestry_excludes2.txt"
rotate_alleles = inpdir + "/lipo_gwas_strict_replication" + "/" + "rotate_alleles.txt"
annotate_vars = inpdir + "/lipo_gwas_strict_replication" + "/" + "gwas_top_hits_annot.txt"
annotate_vars_rsid = inpdir + "/lipo_gwas_strict_replication" + "/" + "gwas_top_hits_annot_rsid.txt"
print("--INFO: Running the job: ")

gwas_hits = [x.split("\t") for x in csvlines2list(snps_file)]
gwas_hits_coords = [[x.split("\t")[0],x.split("\t")[1],x.split("\t")[1]] for x in csvlines2list(snps_file)]
gwas_hits_coords_as_string = ",".join([":".join(x[0:2]) for x in gwas_hits_coords])

vcfs_to_use = [x.split("/home/dgrigoriadis")[1] for x in list(set([x[0] for x in [coordlist2vcf(x, map_df) for x in gwas_hits_coords]]))]

#Extract samples and the gwas hits variants from the agg_v2 vcf files
outvcfs = []
for vcf in vcfs_to_use:
    print(vcf)
    prefix = os.path.basename(vcf.split(".")[0])
    func_vcf = funvcfs+"/"+os.path.basename(vcf.split(".")[0])+"_VEPannot.vcf.gz"
    rjd = project_js + "/" + prefix + ".jobscript"
    time.sleep(0.1)
    scommand = bcftools + ' view '+(' -S '+sample_file+' --force-samples ' if sample_file!="" else "")+\
               '-r '+gwas_hits_coords_as_string+" "+\
               '-Oz -o '+genomic_data+"/"+prefix+"."+project+'.vcf.gz '+vcf+";\n\n"\
               ""+\
               tabix + " " + genomic_data+"/"+prefix+"."+project+'.vcf.gz;\n\n'+\
               ""+\
               "exit"
    #time.sleep(0.1)
    outvcfs.append(genomic_data+"/"+prefix+"."+project+'.vcf.gz')
    #print(scommand)
    lsf_submit(scommand, jobprefix=prefix+"_v2", to_wait_id="",
    wtime="1:00:00", nodes=1, cpu=threads, mem="4gb", cwd=project_js)
    #break

#concat the vcfs all together
scommand2 = bcftools + ' concat -a  '+" ".join(outvcfs)+' '+\
               '-Ov -o '+genomic_data+"/"+project+'.vcf'+";\n\n"+\
               "exit"
time.sleep(0.1)
#print(scommand)
lsf_submit(scommand2, jobprefix=project+"_concat", to_wait_id="",
wtime="1:00:00", nodes=1, cpu=threads, mem="4gb", cwd=project_js)
#break
               ""+\
               "rm "+" ".join(outvcfs)+";\n\n"\

#Run plink commands to create the plink file sets and perform the association analysis
#We need to rename the variants and rotate alleles to match the discovery gwas hits:
scommand3 = plink2 + ' --vcf ' + genomic_data+"/"+project+".vcf --vcf-half-call 'missing'" +\
                     ' --fam ' + fam_file_y5 +\
                     ' --set-missing-var-ids @:#\$r\$a '+\
                     ' --ref-allele force '+rotate_alleles+\
                     ' --make-bed ' +\
                     ' --out ' + genomic_data+"/"+project+"_temp1" +";\n\n"+\
                     ''+\
                     ''+\
            'plink2'+' --bfile '+ genomic_data+"/"+project+"_temp1" +\
                     ' --set-all-var-ids @:#\$r\$a '+\
                     ' --make-bed '+\
                     ' --out '+genomic_data+"/"+project+"_temp2" +";\n\n"+\
            'plink2'+' --bfile '+ genomic_data+"/"+project+"_temp2" +\
                     ' --update-name '+annotate_vars+\
                     ' --extract '+annotate_vars_rsid+\
                     ' --make-bed '+\
                     ' --out '+genomic_data+"/"+project+";\n\n"+\
            'plink2'+' --bfile '+ genomic_data+"/"+project +\
                     ' --keep ' + fam_file_y4 +\
                     ' --remove ' + exclude_ancestry + \
                     ' --glm allow-no-covars cols=chrom,pos,ref,alt1,alt,a1countcc,totallelecc,nobs,gcountcc,beta,orbeta,tz,p'+\
                     ' --out ' + genomic_data+"/"+project+"_y4;\n\n" +\
            plink+' --bfile '+ genomic_data+"/"+project +\
                     ' --keep ' + fam_file_y4 +\
                     ' --remove ' + exclude_ancestry +\
                     ' --logistic --ci 0.95'+\
                     ' --out ' + genomic_data+"/"+project+"_y4;\n\n" +\
            ''+\
            ''+\
            'rm *temp*;\n\n'+\
            'exit'

lsf_submit(scommand3, jobprefix=project+"_plink", to_wait_id="",
wtime="1:00:00", nodes=1, cpu=threads, mem="4gb", cwd=project_js)

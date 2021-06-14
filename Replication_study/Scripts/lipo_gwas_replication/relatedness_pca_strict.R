library(Rlabkey)
library(dplyr)
library(data.table)

#Find samples and controls for GWAS replication study

#Load the samples which have been used in the Discovery Study so we can exclude them
strict_lipo_gids = read.table("../../Input/lipo_gwas_strict_replication/strict_gwas_used_gel_ids.txt",stringsAsFactors = F,header = F)$V1

#Exclude samples
exclude_pids = read.table("../../Input/lipo_gwas_super-strict_replication/exclude_pids.txt",stringsAsFactors = F,header = F)$V1

#Prioritised GEL samples files Y1 to Y5 prioritisation
y1y2_lipo_pids = read.table("../../Input/Y1_Y2_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y1toy5_lipo_pids = read.table("../../Input/Y1-Y5_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y1_lipo = read.table("../../Input/Y1_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y2_lipo = read.table("../../Input/Y2_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y3_lipo = read.table("../../Input/Y3_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y4_lipo = read.table("../../Input/Y4_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y5_lipo = read.table("../../Input/Y5_lipoedema_participant_ids.txt", stringsAsFactors = F, header = F)$V1
y1toy5_lipo = data.frame(rbind(cbind(y1_lipo,"y1"),cbind(y2_lipo,"y2"),cbind(y3_lipo,"y3"),cbind(y4_lipo,"y4"),cbind(y5_lipo,"y5")),stringsAsFactors = F)
colnames(y1toy5_lipo) = c("IID","Y")
y1toy5_lipo_pids = y1toy5_lipo$IID

#Get the labkey entries for all the lipoedema samples prioritised before
y1toy5_labkey = labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v11_2020-12-17", 
  schemaName="lists", 
  queryName="rare_disease_analysis", 
  viewName="", 
  colSelect="participant_id,plate_key,rare_diseases_family_id",
  colFilter=makeFilter(c("participant_id", "IN", paste(y1toy5_lipo_pids,collapse=";"))), 
  containerFilter=NULL, 
  colNameOpt="rname")

#Remove duplicated entries
y1toy5_labkey = y1toy5_labkey[!duplicated(y1toy5_labkey$participant_id),]
y1toy5_all = merge(y1toy5_lipo,y1toy5_labkey,by.x="IID",by.y="participant_id",all.x=T,all.y=F)
write.table(y1toy5_all,"../../Input/lipoedema_ycat_pid_platekey.txt", sep="\t", quote=F, row.names = F, col.names = F)

#Find Lipo samples
##Strict

###Lipoedema samples which have already been used in the discovery study
strict_used_lipocases_temp <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v11_2020-12-17", 
  schemaName="lists", 
  queryName="rare_disease_analysis", 
  viewName="", 
  colSelect="participant_id,plate_key,rare_diseases_family_id,biological_relationship_to_proband,participant_type,participant_phenotypic_sex,participant_ethnic_category,normalised_specific_disease,genome,genome_build,path,tiering,reported_phenotypic_sex,reported_karyotypic_sex,inferred_sex_karyotype,genetic_vs_reported_results", 
  colFilter=makeFilter(c("participant_id", "IN", paste(strict_lipo_gids,collapse=";"))), 
  containerFilter=NULL, 
  colNameOpt="rname"
)

###Lipoedema samples which can be used in the replication study
strict_to_use_lipocases <- labkey.selectRows(
  baseUrl="https://labkey-embassy.gel.zone/labkey", 
  folderPath="/main-programme/main-programme_v11_2020-12-17", 
  schemaName="lists", 
  queryName="rare_disease_analysis", 
  viewName="", 
  colSelect="participant_id,plate_key,rare_diseases_family_id,biological_relationship_to_proband,participant_type,participant_phenotypic_sex,participant_ethnic_category,normalised_specific_disease,genome,genome_build,path,tiering,reported_phenotypic_sex,reported_karyotypic_sex,inferred_sex_karyotype,genetic_vs_reported_results", 
  colFilter=makeFilter(c("rare_diseases_family_id", "NOT_IN", paste(unique(strict_used_lipocases_temp$rare_diseases_family_id), collapse=";")),
                                                                                      c("participant_id", "IN", paste(y1toy5_lipo_pids, collapse=";"))),  
  containerFilter=NULL, 
  colNameOpt="rname"
)


strict_to_use_pk = unique(strict_to_use_lipocases$plate_key)

#Load the masterfile created before
masterfile = read.table("/home/dgrigoriadis/re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/agg_v2/rare_tests/Input/masterfile.txt",sep="\t",stringsAsFactors = F)
rownames(masterfile)=masterfile$IID
colnames(masterfile)=c("IID","FID","Type","Sex","Phenotype")
###Only keep the lipoedema samples defined before.
masterfile = masterfile[(masterfile$IID %in% strict_to_use_pk | masterfile$Phenotype=="Control") & masterfile$Sex=="Female",]
controls = masterfile[masterfile$Phenotype=="Control",]$IID


##Check relatedness and remove related samples
ibd_tab = data.frame(fread("~/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/PCs_relatedness/relatedness/GEL_aggV2_MAF5_mp10_0.0442.kin0"),stringsAsFactors = F)
ibd_tab = ibd_tab[ibd_tab$IID1 %in% unique(c(strict_to_use_pk,super_strict_to_use_pk,controls)) & ibd_tab$IID2 %in% unique(c(strict_to_use_pk,super_strict_to_use_pk,controls)),]
unrelated = data.frame(fread("~/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/PCs_relatedness/relatedness/GEL_aggV2_MAF5_mp10.king.cutoff.unrelated.id"),stringsAsFactors = F)
related = data.frame(fread("~/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/PCs_relatedness/relatedness/GEL_aggV2_MAF5_mp10.king.cutoff.related.id"),stringsAsFactors = F)

related_controls = controls[controls %in% unique(c(ibd_tab$IID1,ibd_tab$IID2))]
ibd_no_con = ibd_tab[!(ibd_tab$IID1 %in% related_controls | ibd_tab$IID2 %in% related_controls),]

#Exclude related cases
exclude_cases = c("XXXXXXXXX-XXX_XXX")

unrelated_masterfile_strict = masterfile[(masterfile$Phenotype=="Control" | masterfile$IID %in% strict_to_use_pk) & !masterfile$IID %in% unique(c(related_controls,exclude_cases)),]
unrelated_masterfile_super_strict = masterfile[(masterfile$Phenotype=="Control" | masterfile$IID %in% super_strict_to_use_pk) & !masterfile$IID %in% unique(c(related_controls,exclude_cases)),]

#Ancestry INFERENCE
ancestries = data.frame(fread("~/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/ancestry/MAF5_superPop_predicted_ancestries.tsv"),stringsAsFactors = F)
eur_ancestries = ancestries[ancestries$eur>=0.9,]

pca_eur = ggplot(pca[pca$IID %in% eur_ancestries$Sample,], aes(x=pc1, y=pc2, label=IID)) + 
  geom_point(size=2, shape=21) +
  labs(x="PC1", y = "PC2") + #geom_text(aes(label=ifelse((abs(PC1)>4*IQR(PC1)|abs(PC2)>4*IQR(PC2)),IID,"")), vjust=-0.6) +
  theme_classic()

pca_eur
ggplotly(pca_eur)

#Create the final samples masterfile
unrelated_masterfile_strict = merge(unrelated_masterfile_strict,y1toy5_all[,c("plate_key","Y")],by.x="IID",by.y="plate_key",all.x=T,all.y=F)
unrelated_masterfile_strict$Y = ifelse(is.na(unrelated_masterfile_strict$Y),"Control",unrelated_masterfile_strict$Y)
###Remove the ethnic outliers
unrelated_masterfile_strict = unrelated_masterfile_strict[(unrelated_masterfile_strict$Phenotype=="Control" & unrelated_masterfile_strict$IID %in% eur_ancestries$Sample) | unrelated_masterfile_strict$Phenotype!="Control",]
unrelated_masterfile_strict$ancestry = "eur"
unrelated_masterfile_strict[unrelated_masterfile_strict$Phenotype!="Control" & !unrelated_masterfile_strict$IID %in% eur_ancestries$Sample,]$ancestry = "exclude"

write.table(unrelated_masterfile_strict,
            "/home/dgrigoriadis/re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/agg_v2/data_requests/Input/lipo_gwas_strict_replication/masterfile.txt",
            quote = F, row.names = F, col.names = F, sep="\t")

write.table(unrelated_masterfile_strict$IID,
            "/home/dgrigoriadis/re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/agg_v2/data_requests/Input/lipo_gwas_strict_replication/masterfile_samples.txt",
            quote = F, row.names = F, col.names = F, sep="\t")


#Create plink fam files
fam_maker <- function(y_masterfile,filename){
  y_masterfile$fam_pheno = ifelse(y_masterfile$Phenotype=="Control",1,2)
  y_masterfile$fam_sex = ifelse(y_masterfile$Sex=="Female",2,1)
  y_masterfile$father_ID = 0
  y_masterfile$mother_ID = 0
  print()
  write.table(y_masterfile[y_masterfile$Y %in% c("Control","y1"),
                           c("FID","IID","father_ID","mother_ID","fam_sex","fam_pheno")],paste0(filename,".y1.fam"),
                           quote = F, row.names = F, col.names = F, sep="\t")
  write.table(y_masterfile[y_masterfile$Y %in% c("Control","y1","y2"),
                           c("FID","IID","father_ID","mother_ID","fam_sex","fam_pheno")],paste0(filename,".y1y2.fam"),
                           quote = F, row.names = F, col.names = F, sep="\t")
  write.table(y_masterfile[y_masterfile$Y %in% c("Control","y1","y2","y3"),
                           c("FID","IID","father_ID","mother_ID","fam_sex","fam_pheno")],paste0(filename,".y1y2y3.fam"),
                            quote = F, row.names = F, col.names = F, sep="\t")
  write.table(y_masterfile[y_masterfile$Y %in% c("Control","y1","y2","y3","y4"),
                           c("FID","IID","father_ID","mother_ID","fam_sex","fam_pheno")],paste0(filename,".y1y2y3y4.fam"),
                           quote = F, row.names = F, col.names = F, sep="\t") 
  write.table(y_masterfile[y_masterfile$Y %in% c("Control","y1","y2","y3","y4","y5"),
                           c("FID","IID","father_ID","mother_ID","fam_sex","fam_pheno")],paste0(filename,".y1y2y3y4y5.fam"),
                           quote = F, row.names = F, col.names = F, sep="\t")  

}

##Strict
fam_maker(unrelated_masterfile_strict,"/home/dgrigoriadis/re_gecip/cardiovascular/postergaard/NGS/Genomes_from_GEL/hg38/agg_v2/data_requests/Input/lipo_gwas_strict_replication/masterfile")

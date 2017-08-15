
### prepare sample list and phenotype file for Folate GxE study in PDB540 and PDB1724
### main tasks: 1) avoid overlap 2) assign one pregnancy to one mother 3) generate pheno file
### 2017 August 10-15. Jonas B.

# THESE INPUT FILES WILL BE NEEDED ON THE GO
# "~/Biostuff/MOBA_PDB1724/interkey_20170721_jb.txt"
# "~/Biostuff/MOBA_EPIDEM_DATA/MoBa_v9/PDB540_SV_INFO_v9.csv"
# "~/Biostuff/MOBA_GENETIC_DATA/MoBa_QCed_5530/H_FINAL/MoBa_5530_EurQC_hg19.fam"
# "~/Biostuff/MOBA_PDB1724/sample_flag_list-20170522.txt"
# "~/Biostuff/MOBA_PDB1724/PCA/pca.eigenvec"
# "~/Biostuff/MOBA_PDB1724/PCA/pca.eigenval"
# "~/Biostuff/MOBA_PDB1724/PCA/PDB540_key_tofakeids_pca.fam"
# "~/Biostuff/MOBA_PDB1724/PCA/igsr_samples.tsv"
# "~/Biostuff/MOBA_PDB1724/Linkage_PDB1724.csv"
# "~/Biostuff/MOBA_PDB1724/Data/PDB1724_MBRN_460_v9.csv"
# "~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q1_v9.csv"
# "~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q2_calculation_v9.csv"
# "~/Biostuff/MOBA_FOLATE_GXE/pdb540_folate_jonas.csv"

# REMAINING TO DO:
# ... think of including Questionnaire version: "VERSJON_KOST_TBL1"
# ... QC shoudl exclude too related individuals. currently this is missing


######################
######################   CHOOSE THE BEST PREGNANCY FOR EACH MOTHER
######################


# load interkey 
interkey = read.table("~/Biostuff/MOBA_PDB1724/interkey_20170721_jb.txt",h=T,stringsAsFactors = F)

####   PROJECT PDB 540
pdb540_inf = read.table("~/Biostuff/MOBA_EPIDEM_DATA/MoBa_v9/PDB540_SV_INFO_v9.csv",h=T,sep=",",stringsAsFactors = F)
pdb540_fam = read.table("~/Biostuff/MOBA_GENETIC_DATA/MoBa_QCed_5530/H_FINAL/MoBa_5530_EurQC_hg19.fam",h=F,stringsAsFactors = F)
pdb540_fam = pdb540_fam[which(substr(pdb540_fam$V2,5,5)=="M"),]
pdb540_fam = pdb540_fam[which(pdb540_fam$V1 %in% pdb540_inf$PREG_ID_540),]  # only legally valid for analysis
mids_pdb540 = pdb540_inf$M_ID_540[which(pdb540_inf$PREG_ID_540 %in% pdb540_fam$V1)] # mother IDs
pregIDs_pdb540 = unique(pdb540_inf$PREG_ID_540[which(pdb540_inf$M_ID_540 %in% mids_pdb540)]) # all pregs of these moms
#pregIDs_pdb540 = pdb540_fam$V1  # it would be a mistake to use this option
reserved_pregIDs_pdb1724 = interkey$PREG_ID_1724[which(interkey$PREG_ID_540 %in% pregIDs_pdb540)] # not to be used in PDB1724
# 388 (correct) vs 312 (not correct) reserved pregIDs

###
###  ONLY SAMPLES WITH GOOD GENOTYPING QUALITY AND CERTAIN IDENTITY
###

pdb1724_flag = read.table("~/Biostuff/MOBA_PDB1724/sample_flag_list-20170522.txt",h=T,stringsAsFactors = F)
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$ROLE=="PARENT"),]
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$select__DUPLICATES=="P"),]
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$phenotypesOK=="TRUE"),] # OK quality, OK Identity

# note (fraction of ethnic outliers in original PCA)
# mean(pdb1724_flag$identify__PCA=="F",na.rm=T)  #  ~ 5%

###
### ONLY SAMPLES WITH HOMOGENEOUS ETHNICITY
###


# load principal components ...
pdb1724_pca = read.table("~/Biostuff/MOBA_PDB1724/PCA/pca.eigenvec",h=F,stringsAsFactors = F)
# ... and their relative importance (weights)
evs = read.table("~/Biostuff/MOBA_PDB1724/PCA/pca.eigenval",h=F,stringsAsFactors = F)

# remove MoBa2008
pca_helper = read.table("~/Biostuff/MOBA_PDB1724/PCA/PDB540_key_tofakeids_pca.fam",h=F,stringsAsFactors = F)
pdb1724_pca = pdb1724_pca[which(!pdb1724_pca$V2 %in% pca_helper$V4),]; rm(pca_helper)
# remove 1000G
pca_helper = read.table("~/Biostuff/MOBA_PDB1724/PCA/igsr_samples.tsv",h=T,sep="\t",stringsAsFactors = F)
pca_helper = pca_helper[,c("Sample.name","Population.code","Superpopulation.code")]
pdb1724_pca = pdb1724_pca[which(!pdb1724_pca$V2 %in% pca_helper$Sample.name),]; rm(pca_helper)

## estimate Euclidean distance of each sample from the "Platonic Norwegian"
evs = evs[1:20,1] # ... evs - introduces weights to PCs based on eigenvalues (importance/variance)
evs = evs / sum(evs)
m = pdb1724_pca
dists = (m$V3 - median(m$V3))^2*evs[1] + (m$V4 - median(m$V4))^2*evs[2] + 
        (m$V5 - median(m$V5))^2*evs[3] + (m$V6 - median(m$V6))^2*evs[4] +
        (m$V7 - median(m$V7))^2*evs[5] + (m$V8 - median(m$V8))^2*evs[6] +
        (m$V9 - median(m$V9))^2*evs[7] + (m$V10 - median(m$V10))^2*evs[8] +
        (m$V11 - median(m$V11))^2*evs[9] + (m$V12 - median(m$V12))^2*evs[10] +
        (m$V13 - median(m$V13))^2*evs[11] + (m$V14 - median(m$V14))^2*evs[12] +
        (m$V15 - median(m$V15))^2*evs[13] + (m$V16 - median(m$V16))^2*evs[14] +
        (m$V17 - median(m$V17))^2*evs[15] + (m$V18 - median(m$V18))^2*evs[16] +
        (m$V19 - median(m$V19))^2*evs[17] + (m$V20 - median(m$V20))^2*evs[18] +
        (m$V21 - median(m$V21))^2*evs[19] + (m$V22 - median(m$V22))^2*evs[20] 
dists = sqrt(dists)

## identify PCA outliers
#mean(dists>0.007)  # suggested fraction of ethnic outliers  (~2%)
pca_outlier_IDs = m$V2[which(dists>0.007)]  # arbitrary  ***
# note: 0.007 threshold is set according to jump in min Pval in Folate association with 20 PCs

# reduce the flag list using this current information about "non-Norwegiousness"
pdb1724_flag = pdb1724_flag[which(!pdb1724_flag$IID %in% pca_outlier_IDs),]


###
### ONLY MOTHERS, not in MoBa2008, good QC, Norwegian
###

pdb1724_link = read.table("~/Biostuff/MOBA_PDB1724/Linkage_PDB1724.csv",h=T,stringsAsFactors = F,sep=",")
pdb1724_link = pdb1724_link[which(!pdb1724_link$PREG_ID_1724 %in% reserved_pregIDs_pdb1724),] # nonoverlap
pdb1724_link = pdb1724_link[which(pdb1724_link$Role=="Mother"),] # only mothers
pdb1724_link = pdb1724_link[which(pdb1724_link$SentrixID_1 %in% pdb1724_flag$IID),] # only good quality
candidate_PregIDs_1724 = pdb1724_link$PREG_ID_1724  # used to reduce folowing datasets

###
### load all relevant datasets
###

# load MFR
mfr_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_MBRN_460_v9.csv",h=T,sep=",",stringsAsFactors = F)
#mean(mfr_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)  
#mean(candidate_PregIDs_1724 %in% mfr_p1724$PREG_ID_1724)  
mfr_p1724 = mfr_p1724[which(mfr_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
mfr_p1724 = mfr_p1724[,c("PREG_ID_1724","FSTART","LEIE","KSNITT","KSNITT_PLANLAGT","VANNAVGANG",
                         "MORS_ALDER","PARITET_5","FLERFODSEL","PLURAL","PLUREK_1",
                         "DODKAT","EKLAMPSI","PREEKL","PREEKL_EKLAMPSI","IVF",
                         "SVLEN_DG","VEKT","KJONN")] #"KSNITT_TIDLIGERE"

# load Q1
Q1_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q1_v9.csv",h=T,sep=",",stringsAsFactors = F)
#mean(Q1_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)
#mean(candidate_PregIDs_1724 %in% Q1_p1724$PREG_ID_1724)
Q1_p1724 = Q1_p1724[which(Q1_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
Q1_p1724 = Q1_p1724[,c("PREG_ID_1724","AA85","AA87")]  # for BMI estimation

# load FFQ
FFQ_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q2_calculation_v9.csv",h=T,sep=",",stringsAsFactors = F)
#mean(FFQ_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)
#mean(candidate_PregIDs_1724 %in% FFQ_p1724$PREG_ID_1724)
FFQ_p1724 = FFQ_p1724[which(FFQ_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
FFQ_p1724 = FFQ_p1724[,c("PREG_ID_1724","KJ","FOLAT")]

# load and add Folate from supplements estimates
fol = read.table("~/Biostuff/MOBA_FOLATE_GXE/pdb540_folate_jonas.csv",h=T,sep=",",stringsAsFactors = F)
fol = fol[which(fol$PREG_ID_540 %in% interkey$PREG_ID_540),]
fol = merge(interkey,fol,by="PREG_ID_540",all=F)

# merge
m2 = merge(mfr_p1724,Q1_p1724,by="PREG_ID_1724",all=T)
m1 = merge(m2,FFQ_p1724,by="PREG_ID_1724",all=T)
m = merge(m1,fol,by="PREG_ID_1724",all.x=T)

#table(FFQ_NA=is.na(m$FOLAT),SUP_NA=is.na(m$folate_food))


######  PREVIEW  #######
# table(m$FSTART,useNA = "a") # 1= spontaneous
# table(m$LEIE,useNA = "a") # 1= normal  (2=breech,3=transverse)
# table(m$KSNITT,useNA = "a") # 1= planned, 2=emergency, 3=unspecified
# table(m$VANNAVGANG,useNA = "a") # 1= 12-24h before delivery,2= >24h before delivery, 3=unspecif
# table(m$DODKAT,useNA = "a") # 0,liveborn,still alive, 11=liveborn, unknown status, 12, liveborn emigrated
# table(m$IVF,useNA = "a") # 1,2,9 = IVF
# table(m$FLERFODSEL,useNA = "a") # 0=singleton,1=plural
# table(m$PLURAL,useNA = "a") # 1=singleton
# table(m$PLUREK_1,useNA = "a") # 1=singleton
# table(m$PREEKL_EKLAMPSI,useNA = "a") # 1=yes
# table(m$AA85,useNA = "a") # 1=yes
# table(m$AA87,useNA = "a") # 1=yes
# table(is.na(m$KJ))
# table(is.na(m$FOLAT))


###  essential data, required
bad00 = which(is.na(m$SVLEN_DG))
bad01 = which(is.na(m$KJ)) #  same as folate is missing
bad_rix = unique(c(bad00,bad01))
m$flag1 = 0
m$flag1[bad_rix] = 1
rm(bad_rix)

### less important, but important (second tier of preference)
bad1 = which(m$FSTART %in% c(2,3))
bad2 = which(m$LEIE %in% c(2,3))
bad3 = which(m$KSNITT %in% c(1))
bad4 = which(m$VANNAVGANG %in% c(1,2,3)) ###   ***
bad5 = which(m$IVF %in% c(1,2,9))
bad6 = which(m$PREEKL_EKLAMPSI == 1)
bad7 = which( (is.na(m$AA85))|(is.na(m$AA87)) )
bad8 = which( (m$AA85<20)|(m$AA85>200) )
bad9 = which( (m$AA87<125)|(m$AA85>200) )
bad10 = which(m$SVLEN_DG > 311)
bad11 = which((m$KJ < 2000)|(m$KJ > 20000))  # hist(m$KJ,breaks=100,col="Grey")
bad12 = which( is.na(m$FLERFODSEL)|(m$FLERFODSEL==1) )
bad13 = which( is.na(m$folate_suppl)&(!is.na(m$FOLAT)) )
bad_rix = unique(c(bad1,bad2,bad3,bad4,bad5,bad6,bad7,bad8,bad9,bad10,bad11,bad12,bad13))
m$flag2 = 0
m$flag2[bad_rix] = 1
rm(bad_rix)

#table(m$flag1,m$flag2)
#table(m$PARITET_5)


### the final prioritisation is  done based on the parity
m$flag3 = m$PARITET_5
m$flag3[which(is.na(m$PARITET_5))] = 5  #  any value higher than 1

#######
####### THE BIG MERGE
#######

inds = pdb1724_link[,c("PREG_ID_1724","SentrixID_1")]
#head(inds); table(table(inds$SentrixID_1)); table(table(inds$PREG_ID_1724))
fams = m[,c("PREG_ID_1724","flag1","flag2","flag3")]  # same mother twice sometimes
#table(table(fams$PREG_ID_1724))
mrg = merge(inds,fams,by="PREG_ID_1724",all=F)  # all.x=T

# assign importance to each flag
mrg$selector = mrg$flag1*100 + mrg$flag2*10 + mrg$flag3*1  ##  giving weights
#table(mrg$selector)

# for each mother choose the pregnancy with the best quality data

library(dplyr)
dd = group_by(mrg,SentrixID_1) %>% arrange(selector) %>% 
        summarise(n=n(),mn=min(selector),mx=max(selector),
                  sl=selector[1],PREG_ID_1724=PREG_ID_1724[1]) %>% ungroup()
# doublecheck
#sum(is.na(dd$mn))
#sum(is.na(dd$mx))

# preview
#sum((dd$n>1)&(dd$mn==dd$mx)) # these were not resolved (n=2)
#sum((dd$n>1)&(dd$mn<dd$mx))  # these were inteligently resolved (n=391)

# example of benefits:
#dd[which(dd$n==3),]
#pregIDs = mrg$PREG_ID_1724[mrg$SentrixID_1=="9970499059_R04C02"]
#m[which(m$PREG_ID_1724 %in% pregIDs),]


dd = dd[,c("SentrixID_1","PREG_ID_1724")]
dd = as.data.frame(dd) # table(table(dd$SentrixID_1)); table(table(dd$PREG_ID_1724))

### final dataset
fin = merge(dd,m,by="PREG_ID_1724",all=F)
# table(fin$flag1,useNA = "a"); table(fin$flag2,useNA = "a")

######################
######################   CLEANING
######################

fin = fin[which(fin$flag1==0),] # gestAge and FFQ answers must be present


bad1 = which(fin$FSTART %in% c(2,3))
bad2 = which(fin$LEIE %in% c(2,3))
bad3 = which(fin$KSNITT %in% c(1))
#bad4 = which(fin$VANNAVGANG %in% c(1,2,3)) # ***
bad5 = which(fin$IVF %in% c(1,2,9))
bad6 = which(fin$PREEKL_EKLAMPSI == 1)
#bad7 = which( (is.na(m$AA85))|(is.na(m$AA87)) )
#bad8 = which( (m$AA85<20)|(m$AA85>200) )
#bad9 = which( (m$AA87<125)|(m$AA85>200) )
#bad10 = which(m$SVLEN_DG > 311)
bad11 = which((fin$KJ < 2500)|(fin$KJ > 20000))  # hist(m$KJ,breaks=100,col="Grey")
bad12 = which( is.na(fin$FLERFODSEL)|(fin$FLERFODSEL==1) )
#bad13 = which( is.na(fin$folate_suppl)&(!is.na(fin$FOLAT)) )
bad_rix = unique(c(bad1,bad2,bad3,bad5,bad6,bad11,bad12)) #bad4,bad7,bad8,bad9,bad10,bad13
fin = fin[-bad_rix,]; rm(bad_rix)


fin = fin[,c("SentrixID_1","PREG_ID_1724","SVLEN_DG","MORS_ALDER","PARITET_5","VANNAVGANG","VEKT","KJONN","KJ","FOLAT","folate_food","folate_suppl")]

# some rearangements
rix_mis = which(is.na(fin$folate_food))
fin$folate_food[rix_mis] = fin$FOLAT[rix_mis]
fin$folate_food_adj = fin$folate_food / fin$KJ * median(fin$KJ) # normalize
fin$folate_total = fin$folate_food_adj*0.6 + fin$folate_suppl # standard formula

# doublecheck
#table(is.na(fin$folate_total))  # there should be some missing values
#hist(fin$KJ,breaks=100,col="grey"); abline(v=2500)
#hist(fin$folate_food,breaks=100,col="grey")
#hist(fin$folate_suppl,breaks=50,col="grey")
#plot(fin$folate_food ~ fin$KJ)
#plot(fin$folate_food_adj ~ fin$KJ)
#summary(lm(fin$folate_food_adj ~ fin$KJ))


####  SAVING AND FILE_NAMING

setwd("~/Dropbox/GIT/FOLATE_GxE/2017/")
report = system("git status",intern = T)
if ( length(grep("modified",report))>0) {
        git_hash = "--???--" # warn user that script does not represent the contents
} else {
git_hash = system("git log --pretty=format:'%h' -n 1",intern = T)  # get recent git commit version
}
time_stamp = substr(Sys.time(),1,10)
file_name = paste("FolGxE_moms_JB_",git_hash,"_n",nrow(fin),"_",time_stamp,".RData",sep="")
save(list="fin",file = paste("~/Biostuff/MOBA_FOLATE_GXE/",file_name,sep=""))



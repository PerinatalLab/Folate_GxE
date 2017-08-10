
### prepare sample list and phenotype file for Folate GxE study in PDB540 and PDB1724
### main tasks: 1) avoid overlap 2) assign one pregnancy to one mother
### 2017 August 10. Jonas B.

# load interkey 
interkey = read.table("~/Biostuff/MOBA_PDB1724/interkey_20170721_jb.txt",h=T,stringsAsFactors = F)


####   PROJECT PDB 540
pdb540_inf = read.table("~/Biostuff/MOBA_EPIDEM_DATA/MoBa_v9/PDB540_SV_INFO_v9.csv",h=T,sep=",",stringsAsFactors = F)
head(pdb540_inf)
pdb540_fam = read.table("~/Biostuff/MOBA_GENETIC_DATA/MoBa_QCed_5530/H_FINAL/MoBa_5530_EurQC_hg19.fam",h=F,stringsAsFactors = F)
pdb540_fam = pdb540_fam[which(substr(pdb540_fam$V2,5,5)=="M"),]
pdb540_fam = pdb540_fam[which(pdb540_fam$V1 %in% pdb540_inf$PREG_ID_540),]  # only legally valid for analysis
head(pdb540_fam)
mids_pdb540 = pdb540_inf$M_ID_540[which(pdb540_inf$PREG_ID_540 %in% pdb540_fam$V1)] # mother IDs
pregIDs_pdb540 = unique(pdb540_inf$PREG_ID_540[which(pdb540_inf$M_ID_540 %in% mids_pdb540)]) # all pregs of these moms
#pregIDs_pdb540 = pdb540_fam$V1  # it would be a mistake to use this option
reserved_pregIDs_pdb1724 = interkey$PREG_ID_1724[which(interkey$PREG_ID_540 %in% pregIDs_pdb540)] # not to be used in PDB1724
# 388 (correct) vs 312 (not correct) reserved pregIDs


pdb1724_flag = read.table("~/Biostuff/MOBA_PDB1724/sample_flag_list-20170522.txt",h=T,stringsAsFactors = F)
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$ROLE=="PARENT"),]
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$phenotypesOK=="TRUE"),]
pdb1724_flag = pdb1724_flag[which(pdb1724_flag$coreOK=="TRUE"),]  # ***
head(pdb1724_flag)

pdb1724_link = read.table("~/Biostuff/MOBA_PDB1724/Linkage_PDB1724.csv",h=T,stringsAsFactors = F,sep=",")
head(pdb1724_link); dim(pdb1724_link)
pdb1724_link = pdb1724_link[which(!pdb1724_link$PREG_ID_1724 %in% reserved_pregIDs_pdb1724),] # reduce
pdb1724_link = pdb1724_link[which(pdb1724_link$Role=="Mother"),]
pdb1724_link = pdb1724_link[which(pdb1724_link$SentrixID_1 %in% pdb1724_flag$IID),] # only good quality
head(pdb1724_link); dim(pdb1724_link)
candidate_PregIDs_1724 = pdb1724_link$PREG_ID_1724  # not essential, used to reduce folowing datasets



##### which pregnancies are valid for analysis

# load MFR
mfr_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_MBRN_460_v9.csv",h=T,sep=",",stringsAsFactors = F)
mean(mfr_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)  # 83.30%
mean(candidate_PregIDs_1724 %in% mfr_p1724$PREG_ID_1724)  # 99.98%
mfr_p1724 = mfr_p1724[which(mfr_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
mfr_p1724 = mfr_p1724[,c("PREG_ID_1724","FSTART","LEIE","KSNITT","KSNITT_PLANLAGT","VANNAVGANG",
                         "MORS_ALDER","PARITET_5","FLERFODSEL","PLURAL","PLUREK_1",
                         "DODKAT","EKLAMPSI","PREEKL","PREEKL_EKLAMPSI","IVF",
                         "SVLEN_DG","VEKT","KJONN")] #"KSNITT_TIDLIGERE"

# load Q1
Q1_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q1_v9.csv",h=T,sep=",",stringsAsFactors = F)
mean(Q1_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)  # 83.30%
mean(candidate_PregIDs_1724 %in% Q1_p1724$PREG_ID_1724)  # 99.99%
Q1_p1724 = Q1_p1724[which(Q1_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
Q1_p1724 = Q1_p1724[,c("PREG_ID_1724","AA85","AA87")]

# load FFQ
FFQ_p1724 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q2_calculation_v9.csv",h=T,sep=",",stringsAsFactors = F)
mean(FFQ_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724)  # 83.26%
mean(candidate_PregIDs_1724 %in% FFQ_p1724$PREG_ID_1724)  # 96.44%
FFQ_p1724 = FFQ_p1724[which(FFQ_p1724$PREG_ID_1724 %in% candidate_PregIDs_1724),]
FFQ_p1724 = FFQ_p1724[,c("PREG_ID_1724","KJ","FOLAT")]

m1 = merge(mfr_p1724,Q1_p1724,by="PREG_ID_1724",all=T)
m = merge(m1,FFQ_p1724,by="PREG_ID_1724",all=T)
head(m)



table(m$FSTART,useNA = "a") # 1= spontaneous
table(m$LEIE,useNA = "a") # 1= normal  (2=breech,3=transverse)
table(m$KSNITT,useNA = "a") # 1= planned, 2=emergency, 3=unspecified
table(m$VANNAVGANG,useNA = "a") # 1= 12-24h before delivery,2= >24h before delivery, 3=unspecif
table(m$DODKAT,useNA = "a") # 0,liveborn,still alive, 11=liveborn, unknown status, 12, liveborn emigrated
table(m$IVF,useNA = "a") # 1,2,9 = IVF
table(m$FLERFODSEL,useNA = "a") # 0=singleton,1=plural
table(m$PLURAL,useNA = "a") # 1=singleton
table(m$PLUREK_1,useNA = "a") # 1=singleton
table(m$PREEKL_EKLAMPSI,useNA = "a") # 1=yes
table(m$AA85,useNA = "a") # 1=yes
table(m$AA87,useNA = "a") # 1=yes
table(is.na(m$KJ))
table(is.na(m$FOLAT))


###  essential
table(is.na(m$SVLEN_DG))
table(is.na(m$FOLAT))
bad00 = which(is.na(m$SVLEN_DG))
bad01 = which(is.na(m$KJ)) #  same as folate is missing
bad_rix = unique(c(bad00,bad01))
m$flag1 = 0
m$flag1[bad_rix] = 1
rm(bad_rix)

### less important
bad1 = which(m$FSTART %in% c(2,3))
bad2 = which(m$LEIE %in% c(2,3))
bad3 = which(m$KSNITT %in% c(1))
bad4 = which(m$VANNAVGANG %in% c(1,2,3))
bad5 = which(m$IVF %in% c(1,2,9))
bad6 = which(m$PREEKL_EKLAMPSI == 1)
bad7 = which( (is.na(m$AA85))|(is.na(m$AA87)) )
bad8 = which( (m$AA85<20)|(m$AA85>200) )
bad9 = which( (m$AA87<125)|(m$AA85>200) )
bad10 = which(m$SVLEN_DG > 311)
bad11 = which(m$KJ > 25000)  # hist(m$KJ,breaks=100,col="Grey")
bad12 = which( is.na(m$FLERFODSEL)|(m$FLERFODSEL==1) )
bad_rix = unique(c(bad1,bad2,bad3,bad4,bad5,bad6,bad7,bad8,bad9,bad10,bad11,bad12))
m$flag2 = 0
m$flag2[bad_rix] = 1
rm(bad_rix)

table(m$flag1,m$flag2)
table(m$PARITET_5)



# big merge

inds = pdb1724_link[,c("PREG_ID_1724","SentrixID_1")]
#head(inds)
#table(table(inds$SentrixID_1))
#table(table(inds$PREG_ID_1724))
fams = m[,c("PREG_ID_1724","PARITET_5","flag1","flag2")]
#table(table(fams$PREG_ID_1724))
mrg = merge(inds,fams,by="PREG_ID_1724",all.x=T)
head(mrg)

mrg$selector = mrg$flag1*100 + mrg$flag2*10 + mrg$PARITET_5*1

dd = group_by(mrg,SentrixID_1) %>% arrange(selector) %>% summarise(n=n(),mn=min(selector),sl=selector[1],prID=PREG_ID_1724[1]) %>% ungroup()

table(dd$mn==dd$sl)

### prefered



table(table(pdb1724_link$SentrixID_1))



head(pdb1724_link)
sum(duplicated(pdb1724_link$SentrixID_1))
attach variable indicating which pregnancy to keep



q2 = read.table("~/Biostuff/MOBA_EPIDEM_DATA/MoBa v6/Data/Q2_beregning_PDB540_v6.csv",h=T,stringsAsFactors = F,sep=",")
dim(q2)
mean(pdb540_fam$V1 %in% q2$PREG_ID_540)
mean(pdb540_fam$V1 %in% fol$PREG_ID_540)
mean(q2$PREG_ID_540 %in% fol$PREG_ID_540 )



head(pdb540_fam)
head(interkey)


fol = read.table("~/Biostuff/MOBA_FOLATE_GXE/pdb540_folate_jonas.csv",sep=",",h=T,stringsAsFactors = F)
head(fol)



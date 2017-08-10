library(boot)
library(dplyr)
library(data.table)

#### This code implements semi-parametric bootstraping for the calculation of p-values in any GxE analysis, 
# using linear regression, and when assumptions are thought not to be met. Can also be used for log regression

# According to: Bůžková P, et al., Annals of Human Genetics (2011)

# Steps (G= SNP, E= environment, Y= dependent var, Y^: fitted dependent var):

# 1. Calculate the null hypothesis:
  # Y= B0 + B1xG + B2xE
# 2. Save fitted values Y^.
# 3. Resample Y^ with replacement (bootstrap). 1000 bootstraps for >1 SNP. 10000 for 1 SNP.
# 4. Use resampled Y^ in:
  # Y^= B0 + B1xG + B2xE + B3xGxE
# 5. Save model f statistic
# 6. Compute original f statistic for the interaction model from the original data by fitting:
 # Y= B0 + B1xG + B2xE + B3xGxE
# 7. Compute the p-value by comparing test statistic in step 6 to the distribution in step 5.

# We need a data frame with the following items for each women:
# G for each SNP id
# Dependent variable 
# E 

# Open data

mf= readRDS(file="XXX")

# Save column names in vector according to col number

myvars= colnames(mf)[2:15] # modify 2:15 according to SNP columns number

# for each SNP, we compute the fitted response for the null hypothesis (B3xGxE= 0) model, 
# resample Y^, and use in a model with the interaction term. Save model test statistic (1000 models)
# Steps 1 to 5

outs = NULL
for(i in myvars){
  f <- paste("Y", "~", i,"*E") # change Y= dependent variable name, change E= environment var name
  noint= lm(as.formula(f), data=mf,na.action= na.exclude)
for(c in 1:1000){
  mf$y= sample(fitted(noint), replace = TRUE)
    f1 <- paste("y", "~", i,"*E") # change only E= environment var name 
  mod= lm(as.formula(f1), data=mf,na.action= na.exclude)
  out = summary(mod)$fstatistic[1]
  out$id = paste(i)
    outs = bind_rows(outs, out)
}
}

# remove for clearness

rm(c,f,f1,mod,noint,out,i)

# Compute original f statistic for the interaction model from the original data. Step 6

inter_pval=lapply( mf[,c(2:15)], function(x) summary(lm(mf$SVLEN_DG ~ x))$fstatistic[1]) # modify 2:15 
# according to SNP columns number

# Merge test statistic from original GxE model, with that of the bootstraped models 

interactions= NULL
interactions=cbind(interactions, inter_pval)
interactions=cbind(id = rownames(interactions), interactions)
interactions=as.data.frame(interactions)
names(outs)[names(outs) == 'value'] <- 'nullhip_pval'

outs=merge(outs, interactions, by="id")

# Compare original test statistiscs with bootstraped test statistics = p-value

df= outs %>%
  group_by(id, inter_pval<nullhip_pval) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

# Clear dataframe 

df= df[grep("TRUE", df$`inter_pval < nullhip_pval`), ]
df=df[,c("id", "freq")]

rm(interactions, outs, inter_pval, myvars)

################### Reproducible example with real data (fake hypothesis)######################

# Interactions between B12 SNPs and maternal age on GA at delivery

# Open data

mf= readRDS(file="/HARVEST/MFR + Harv B12 Mothers.Rda")

names(mf)

#  dataset with bimodal distribution outcome 
mf=mf[mf$SVLEN_DG<253 | mf$SVLEN_DG>273,]

# for each SNP, we compute the fitted response for the null hypothesis (B3xGxE= 0) model, 
# resample Y^, and use in a model with the interaction term. Save model test statistic (1000 models)
# Steps 1 to 5

myvars= colnames(mf)[2:8]

outs = NULL
for(i in myvars){
  f <- paste("SVLEN_DG", "~", i,"*MORS_ALDER") 
  noint= lm(as.formula(f), data=mf,na.action= na.exclude)
  for(c in 1:1000){
    mf$y= sample(fitted(noint), replace = TRUE)
    f1 <- paste("y", "~", i,"*MORS_ALDER") 
    mod= lm(as.formula(f1), data=mf,na.action= na.exclude)
    out = summary(mod)$fstatistic[1]
    out$id = paste(i)
    outs = bind_rows(outs, out)
  }
}

rm(c,f,f1,mod,noint,out,i)

inter_pval=lapply( mf[,c(2:8)], function(x) summary(lm(mf$SVLEN_DG ~ x))$fstatistic[1])


interactions= NULL
interactions=cbind(interactions, inter_pval)
interactions=cbind(id = rownames(interactions), interactions)
interactions=as.data.frame(interactions)
names(outs)[names(outs) == 'value'] <- 'nullhip_pval'

outs=merge(outs, interactions, by="id")

df= outs %>%
  group_by(id, inter_pval<nullhip_pval) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))


df= df[grep("TRUE", df$`inter_pval < nullhip_pval`), ]
df=df[,c("id", "freq")]

rm(interactions, outs, inter_pval, myvars)

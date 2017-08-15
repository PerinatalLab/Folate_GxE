library(dplyr)
library(ggplot2)

#### This code implements semi-parametric bootstraping for the calculation of p-values in any GxE analysis, 
# using linear regression, and when assumptions are thought not to be met. 
# In our case, bimodal outcome.
# 2017 August 14. Pol SN

# According to: Bůžková P, et al., Annals of Human Genetics (2011), with slight modification (use resampled
# residuals instead of fitted values). 
# More in Bůžková P, Epidemiologic methods, 2016. Resampled residuals requires G-E independence.

# Steps (G= SNP, E= environment, Y= dependent var, Y*: Y fixing null hypothesis model):

# 1. Calculate the null hypothesis:
# Y= B0 + B1xG + B2xE
# 2. Calculate interaction and save B3 t-statistic 
# Y= B0 + B1xG + B2xE + B3xGxE
# 3. Resample residuals from model in 1 with replacement (bootstrap). 
# n times of resample?
# 4. Use resampled Y* in:
# Y*= B0 + B1xG + B2xE + B3xGxE
# 5. Save B3 t-statistic
# 6. Compute the p-value by comparing the original B3 t-statistic in step 2 to the t-statistics obtained in
# step 5.

# We need a data frame with the following items for each women:
# G for each SNP id
# Dependent variable 
# E 

# Open data

mf= readRDS(file="XXX")

# Save column names in vector according to col number

myvars= colnames(mf)[2:15] # modify 2:15 according to SNP columns number

# Compute the B3 t-statistic for model Y= B0 + B1xG + B2xE + B3xGxE 
inter_tstat=apply(mf[,c(2:15)],2, #2:15 define interval of snps to be tested 
                  function(x) summary(lm(mf$Y ~ x*E))$coefficients[4,3]) # Y= outcome; E= environmental var
df=as.data.frame(inter_tstat)
df=cbind(snpid = rownames(df), df)

# for each SNP, we compute the residuals for the null hypothesis (B3xGxE= 0) model, 
# resample the residuals, and use as outcome in a linear model with the interaction term. 
# Add original B3 t-statistic previously obtained for each corresponding snp and for each resample
# calculate bootstraped p-value = relative freq of original t-statistic> bootstraped t-statistic

outs = NULL
out= NULL
outs2= NULL
for(i in myvars){
  f <- paste("Y", "~", i,"+E") # Y= outcome; E= environmental var
  noint= lm(as.formula(f), data=mf,na.action= na.exclude)
  for(c in 1:10000){ # resampling n? 1000, 10000?
    mf$res= sample(resid(noint),replace = TRUE)
    f1 <- paste("res", "~", i,"*E") # E= environmental var
    mod= lm(as.formula(f1), data=mf,na.action= na.exclude) # fit null hypothesis model 
    out$snpid = paste(i)
    out$tvalue = summary(mod)$coefficients[4,3]  # retreive t-statistic //should we use p-values instead?
    a1=merge(out, df, by="snpid") # add the t-statistic for the GxE term from original data
    a1=a1[,c("tvalue", "inter_tstat")]
    outs= bind_rows(outs, a1)
    rm(a1)
  }
  a3= outs %>%
    summarise(pvalue = mean((abs(as.numeric(inter_tstat)) < abs(as.numeric(tvalue))))) #compare the 2 t-stats
  a3$snpid = paste(i)
  outs2= bind_rows(outs2,a3) # This is the end file we will be interested in
  outs= NULL
}

# remove for clearness

rm(a3, c,f,f1,i,inter_tstat,mod, noint, out, outs, myvars)

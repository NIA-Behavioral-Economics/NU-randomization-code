# Code to think about covariate-balanced randomization
# Programmer: Lucia Petito
# Date: February 19, 2020
# Reference: https://journals.sagepub.com/doi/full/10.1177/0962280211416038
# Demonstrates t-test is robust, even with non-normally distributed continuous variables

# Functions to describe balance -------------------------------------------

cat.balance <- function(df, rand, varname, Nsim, nthresh=2){
  # Inputs
  # df - dataset
  # rand - matrix with each column corresponding to a randomization scheme
  # varname - name of categorical variable from 'df'
  # Nsim - number of simulations
  # nthresh - absolute difference between categories. Default of 2
  r <- rep(NA, Nsim)
  for (i in 1:Nsim){
    ran <- rand[,i]
    temp <- data.frame(subset(df, select=varname), rand = ran)
    varxrand <- table(temp)
    r[i] <- max(abs(varxrand[,1] - varxrand[,2])) <= nthresh
  }
  return(r)
  # Returns a vector describing which simulation scenarios satisfy criteria
}

cont.balance <- function(df, rand, varname, Nsim, thresh=0.5){
  # Inputs
  # df - dataset
  # rand - matrix with each column corresponding to a randomization scheme
  # varname - name of continuous variable from 'df'
  # Nsim - number of simulations
  # thresh - maximum t-statistic from t-test. Default of 0.5
  r <- rep(NA, Nsim)
  for (i in 1:Nsim){
    ran <- rand[,i]
    temp <- data.frame(subset(df, select=varname), rand = ran)
    tstat <- t.test(eval(parse(text=paste0(varname, "~ rand"))), data=temp)$statistic
    r[i] <- abs(tstat) <= thresh
  }
  return(r)
  # Returns a vector describing which simulation scenarios satisfy threshold criteria
}

# Read in data; do basic cleaning -----------------------------------------

# Set to location of data
setwd("R:/Medicine/GIM/Persell/BEAGLE/Common/R33 documents/Stats/data/")

df <- read.csv("simulated_data_example.csv", header=T)

# Create potential randomizations -----------------------------------------
set.seed(1234)
Nsim <- 1000
n <- nrow(df)
toassign <- factor(rep(c("A", "B"), length.out=n))

# Each column  in 'possiblerand' describes a randomization scheme
possiblerand <- replicate(n=Nsim, 
                          expr=sample(toassign, size=n, replace=F), 
                          simplify = 'matrix')


# Check balance in simulated randomization schemes ------------------------

tau <- 0.5 # Threshold for t-tests

# Number of physicians - continuous
phys.balance <- cont.balance(df=df, rand=possiblerand, varname='n_prv', Nsim=Nsim,  thresh=tau)
table(phys.balance)

# UTI treatment - continuous
uaucdeno.balance <- cont.balance(df=df, rand=possiblerand, varname='deno_uauc', Nsim=Nsim,  thresh=tau)
table(uaucdeno.balance)
uaucpct.balance <- cont.balance(df=df, rand=possiblerand, varname='pct_uaucOrdered_woDx', Nsim=Nsim,  thresh=tau)
table(uaucpct.balance)

# Diabetes Mellitus variables - continuous
dmdeno.balance <- cont.balance(df=df, rand=possiblerand, varname='deno_dm', Nsim=Nsim,  thresh=tau)
table(dmdeno.balance)
dmpct.balance <- cont.balance(df=df, rand=possiblerand, varname='pct_overtreated_dm', Nsim=Nsim,  thresh=tau)
table(dmpct.balance)

# PSA variables - continuous
psadeno.balance <- cont.balance(df=df, rand=possiblerand, varname='psa_deno', Nsim=Nsim,  thresh=tau)
table(psadeno.balance)
psapct.balance <- cont.balance(df=df, rand=possiblerand, varname='pct_PSA_ordered', Nsim=Nsim,  thresh=tau)
table(psapct.balance)

# Region - categorical
region.balance <- cat.balance(df=df, rand=possiblerand, varname='Region', Nsim=Nsim, nthresh=2)
table(region.balance)

# Compile balance results into one matrix
compiledtable <- cbind(region.balance,
                       psadeno.balance,
                       psapct.balance,
                       dmdeno.balance,
                       dmpct.balance,
                       uaucdeno.balance,
                       uaucpct.balance,
                       phys.balance)
# See which randomization schemes worked
collapsedtable <- apply(compiledtable, 1, sum)
table(collapsedtable)
sum(collapsedtable==ncol(compiledtable)) 


# Describe a simulation setting that worked --------------------------------------

library(tableone)
# Indices of randomization schemes that worked
goodRandomizations <- which(collapsedtable==ncol(compiledtable)) 

exampledf <- df
exampledf$rand <- possiblerand[,goodRandomizations[1]] # Change '1' to view other possibilities

vars <- c("n_prv", "deno_uauc", "pct_uaucOrdered_woDx", "deno_dm", "pct_overtreated_dm", "psa_deno", "pct_PSA_ordered")

tableOneRand <- CreateContTable(vars = vars, strata = c('rand'), data = exampledf)
print(tableOneRand)


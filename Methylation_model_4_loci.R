#########################################################
### 	 	 RUN BAYESIAN MODEL ON FOUR LOCI 		  ###
###                    			                      ###
#########################################################

### --- Load Libraries --- ###
library(ggplot2)
library(rstan)
library(ggridges)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("C:/Users/c_cro/Documents/PhD/Epigenetics/Methylation/results")
source("plotPost.R")
options(mc.cores = parallel::detectCores())

### --- Load Data File --- ###
orca_all = read.table("orcaData.csv", header = TRUE, sep = ",")

### --- Remove extra loci and sites located downstream from the TSS --- ###
orca4 = orca_all[which(orca_all$locus == "CRF" | orca_all$locus == "ACTB" | orca_all$locus == "BDNF" | orca_all$locus == "NP"),]
orca = orca4[-which(orca4$site == "BDNF_106" | orca4$site == "BDNF_145" | orca4$site == "BDNF_154" | orca4$site == "BDNF_167" | orca4$site == "BDNF_172" | orca4$site == "BDNF_192" | orca4$site == "BDNF_208" | orca4$site == "BDNF_213" |
                      orca4$site == "NP_196" | orca4$site == "CRF_135" | orca4$site == "CRF_148" | orca4$site == "CRF_236" | orca4$site == "CRF_242" | orca4$site == "CRF_257" | orca4$site == "CRF_260" | orca4$site == "CRF_265" | 
                      orca4$site == "CRF_270" | orca4$site == "CRF_286" | orca4$site == "CRF_305" | orca4$site == "CRF_313" | orca4$site == "CRF_319" | orca4$site == "CRF_321"  ),]

orca$locus = droplevels(orca$locus)
orca$site = droplevels(orca$site)
orca$ind = droplevels(orca$ind)


########################### Prepare the data for Stan ###########################  

### --- Standardize AGE --- ###
age     = orca$age
ageMean = mean(age)
ageSD   = sd(age)
zage    = (age - ageMean) / ageSD

### --- Standardize Percent Methylation (C_Ratios) --- ###
permeth     = orca$C_Ratio
permethMean = mean(permeth)
permethSD   = sd(permeth)
zpermeth    = (permeth - permethMean) / permethSD
N = length(permeth)

### --- Categorical Sex Data --- ###
sex = as.numeric(orca$Sex)
sexNames = levels(orca$Sex)
nsexLevels = length(unique(orca$Sex))

### --- Categorical Population Data --- ###
pop = as.numeric(orca$Population)
popNames = levels(orca$Population)
npopLevels = length(unique(orca$Population))

### --- Categorical Site Data --- ###
site = as.numeric(as.factor(orca$site))
siteNames = levels(orca$site)
nsiteLevels = length(unique(orca$site))

### --- Categorical Locus Data --- ###
locus = as.numeric(orca$locus)
locusNames = levels(orca$locus)
nlocusLevels = length(unique(orca$locus))

### --- Categorical Locus Data --- ###
ind = as.numeric(orca$ID)
indNames = levels(orca$ID)
nindLevels = length(unique(orca$ID))


########################### Create a data list for Stan ###########################

dataList = list(
  permeth = zpermeth,
  pop = pop,
  npopLevels = npopLevels,
  sex = sex,
  nsexLevels = nsexLevels,
  age = zage,
  N = N,
  site = site,
  nsiteLevels = nsiteLevels,
  ind = ind,
  nindLevels = nindLevels
)


############################ Define the model & write a string for Stan ###########################

modelstring = "
  data {
    int N;                  // Sample size
    int npopLevels;         // Number of populations
    int nsiteLevels;        // Number of sites
  	int nsexLevels;         // Number of sexes
    int nindLevels;         // Number of individuals
    vector[N] permeth;      // Vector of percent methylation
	  vector[N] age;     	    // Vector of age
    int pop[N];             // Vector of indicators of which population each sample is from
    int site[N];            // Vector of indicators of which site each sample is from
    int sex[N];             // Vector of indicators of which sex each sample is from
    int ind[N];             // Vector of indicators of which individual each sample is from
} 

  parameters {
    real bage;                    // Coefficient for effect of age
    real bsex[nsexLevels];        // Coefficients for effects of being in each sex
    real bpop[npopLevels];        // Coefficients for effects of being in each population     
    real bsite[nsiteLevels];      // Coefficients for effects of each site
    real bind[nindLevels];        // Coefficients for effects of each individual
    real bpopsite[npopLevels, nsiteLevels];        // Coefficients for interaction between population and site
    real sexMean;                 // Mean effect of sex across all samples
    real <lower=0> sexMeanSD;     // sd of the distribution of the effects across all sexes
    real popMean;                 // Mean effect of population across all samples
    real <lower=0> popMeanSD;     // sd of the distribution of the effects across all populations
  	real <lower=0> sigma[npopLevels];               // Coefficients for overall sd
    real indMean;                 // Mean effect of individual across all samples
    real <lower=0> indMeanSD;     // sd of the distribution of the effects across all individuals
    real siteMean;                 // Mean effect of site across all samples
    real <lower=0> siteMeanSD;     // sd of the distribution of the effects across all sites

}

  model {
    // Definitions
    vector[N] mu;

    // Likelihood
    for (i in 1:N) {
        mu[i] = bage*age[i] + bpop[pop[i]] + bsex[sex[i]] + bind[ind[i]] + bsite[site[i]] + bpopsite[pop[i], site[i]];
        permeth[i] ~ normal(mu[i], sigma[pop[i]]);  
    }

    // Priors
    bage ~ cauchy(0, 1);    // Using a cauchy distribution because different life histories could 
                            // create a lot of outliers in methylation patterns
	
  	for (q in 1:npopLevels) { 
      sigma[q] ~ cauchy(0, 1);
    }
  	
  	for (m in 1:nsexLevels) {
       bsex[m] ~ normal(sexMean, sexMeanSD);
    }
	
  	for (n in 1:nindLevels) {
        bind[n] ~ normal(indMean, indMeanSD);
	  }
    
  	for (n in 1:npopLevels) {
        bpop[n] ~ normal(popMean, popMeanSD);
  	}
  	
  	for (z in 1:nsiteLevels) {
        bsite[z] ~ normal(siteMean, siteMeanSD);
  	}
    
    for (q in 1:npopLevels){
      for (r in 1:nsiteLevels){
        bpopsite[q,r] ~ normal(0, 1);
      }
    }
	
    
    // Hyperpriors
    popMean     ~ normal(0, 1);
    popMeanSD   ~ cauchy(0, 1);
	  sexMean     ~ normal(0, 1);
    sexMeanSD   ~ cauchy(0, 1);
    indMean     ~ normal(0, 1);
    indMeanSD   ~ cauchy(0, 1);
    siteMean    ~ normal(0, 1);
    siteMeanSD  ~ cauchy(0, 1);
   }

  generated quantities {
    // Posterior Predictive Variable Definitions
    vector[N] mu_pred;
    vector[N] permeth_pred;

    // For Posterior Predictive Calculations
    for (i in 1:N) {
      mu_pred[i] = bage*age[i] + bpop[pop[i]] + bsex[sex[i]] + bind[ind[i]] + bsite[site[i]] + bpopsite[pop[i], site[i]];
      permeth_pred[i] = normal_rng(mu_pred[i], sigma[pop[i]]);
    }
  }
"
writeLines(modelstring, con = "orca_site_model.stan")


################################################
###              RUN STAN                    ###
################################################
stanFit <- stan(file = "orca_site_model.stan", 
                data = dataList, 
                pars = c("bage", "bpop", "bsex", "bind", "bsite", "bpopsite", "sigma", "permeth_pred"),
                warmup = 2000,
                iter = 12000, 
                chains = 3)


############################ Check MCMC Performance ############################

print(stanFit)  
stan_trace(stanFit, pars = c("bage", "bsex", "bpop"), inc_warmup = TRUE)
stan_trace(stanFit, pars = "bind", inc_warmup = TRUE)
stan_trace(stanFit, pars = "bpopsite", inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)


############################ Plot Results With Stan Functions ############################

stan_plot(stanFit, par = c("bage", "bsex", "bpop"))
stan_plot(stanFit, par = "bind")
stan_plot(stanFit, par = "bsite")
stan_plot(stanFit, par = "sigma")
stan_plot(stanFit, par = "bpopsite")

s <- stan_plot(stanFit, par = c("bage", "bsex", "bpop"))
s + scale_y_discrete(labels=c("Age","Female","Male","Northern Resident","Southern Resident"))
s

############################ Posterior Predictive Check ############################

### --- Extract the predictions --- ###
mcmcChains = as.data.frame(stanFit)
write.csv(mcmcChains,'mcmcChains.csv')

chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("permeth_pred[", i, "]", sep = "")]
}

### --- Mean expected value for record --- ###
ypredMean = apply(zypred, 2, mean)

### --- Upper and lower expected 95% HDI for each visit --- ###
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)


############################ Model Fit Plots ############################         

### --- Plot mean predicted values --- ###
par(mfrow = c(1, 1))
record = 1:500

subsample = sample(1:N, size = 500, replace=FALSE)
dotchart(ypredMean[subsample], xlim = c(-6, 6), xlab = "Standardized Percent Methylation", ylab = "Sample")

### --- Add HDIs --- ###
segments(x0 = ypredLow[subsample], y0 = record, x1 = ypredHigh[subsample], y1 = record)
### --- Add observed values --- ###
points(x = zpermeth[subsample], y = record, pch = 16, col = rgb(0, 0.8, 1, 0.6))


################################################################
###               PLOTS TO ASSESS EFFECTS                    ###
################################################################

############################ Plot Results With Stan Functions ############################
############################         For Site Effects         ############################

actb<-c(which(startsWith(siteNames,"ACTB")))
bdnf<-c(which(startsWith(siteNames,"BDNF")))
crf<-c(which(startsWith(siteNames,"CRF")))
np<-c(which(startsWith(siteNames,"NP")))

actb_sites<-c()
for (i in actb){
  actb_sites<-append(actb_sites,paste("bsite[",i,"]", sep=""))
}

bdnf_sites<-c()
for (i in bdnf){
  bdnf_sites<-append(bdnf_sites,paste("bsite[",i,"]", sep=""))
}

crf_sites<-c()
for (i in crf){
  crf_sites<-append(crf_sites,paste("bsite[",i,"]", sep=""))
}

np_sites<-c()
for (i in np){
  np_sites<-append(np_sites,paste("bsite[",i,"]", sep=""))
}

stan_plot(stanFit, par = c(actb_sites))
stan_plot(stanFit, par = c(bdnf_sites))
stan_plot(stanFit, par = c(crf_sites))
stan_plot(stanFit, par = c(np_sites))


############################ Plot Results With Stan Functions ############################
############################        For Pop*Site Effects      ############################


actb_ps<-c()
for (i in actb){
  actb_ps<-append(actb_ps,paste("bpopsite[1,",i,"]", sep=""))
  actb_ps<-append(actb_ps,paste("bpopsite[2,",i,"]", sep=""))
}

bdnf_ps<-c()
for (i in bdnf){
  bdnf_ps<-append(bdnf_ps,paste("bpopsite[1,",i,"]", sep=""))
  bdnf_ps<-append(bdnf_ps,paste("bpopsite[2,",i,"]", sep=""))
}

crf_ps<-c()
for (i in crf){
  crf_ps<-append(crf_ps,paste("bpopsite[1,",i,"]", sep=""))
  crf_ps<-append(crf_ps,paste("bpopsite[2,",i,"]", sep=""))
}

np_ps<-c()
for (i in np){
  np_ps<-append(np_ps,paste("bpopsite[1,",i,"]", sep=""))
  np_ps<-append(np_ps,paste("bpopsite[2,",i,"]", sep=""))
}

stan_plot(stanFit, par = c(actb_ps))
stan_plot(stanFit, par = c(bdnf_ps))
stan_plot(stanFit, par = c(crf_ps))
stan_plot(stanFit, par = c(np_ps))

plot(stanFit, par = c(actb_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(bdnf_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(crf_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(np_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")




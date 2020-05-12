#-----------------------#
# Read the data into R  #
#-----------------------#
library(ggplot2)
library(rstan)
library(ggridges)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("C:/Users/c_cro/Documents/PhD/Epigenetics/Methylation/results")
source("plotPost.R")
options(mc.cores = parallel::detectCores())
orca = read.table("orca_4_loci.csv", header = TRUE, sep = ",")

#------------------------------#
# Prepare the data for Stan    #
#------------------------------#
#--- Standardize AGE ---#
age     = orca$age
ageMean = mean(age)
ageSD   = sd(age)
zage    = (age - ageMean) / ageSD

#--- Standardize METHYLATION C_Ratios ---#
permeth     = orca$C_Ratio
permethMean = mean(permeth)
permethSD   = sd(permeth)
zpermeth    = (permeth - permethMean) / permethSD
N = length(permeth)

#---Categorical Sex Data  --#
sex = as.numeric(orca$Sex)
sexNames = levels(orca$Sex)
nsexLevels = length(unique(orca$Sex))

#--- Categorical Population Data ---#
pop = as.numeric(orca$Population)
popNames = levels(orca$Population)
npopLevels = length(unique(orca$Population))

#--- Categorical Site Data ---#
site = as.numeric(as.factor(orca$site))
siteNames = levels(orca$site)
nsiteLevels = length(unique(orca$site))

#--- Categorical Locus Data ---#
locus = as.numeric(orca$locus)
locusNames = levels(orca$locus)
nlocusLevels = length(unique(orca$locus))

#--- Categorical Locus Data ---#
ind = as.numeric(orca$ID)
indNames = levels(orca$ID)
nindLevels = length(unique(orca$ID))


#-----------------------------#
# Create a data list for STAN #
#-----------------------------#
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


#-----------------------------#
# Define the model and write  #
# a string for Stan           #
#-----------------------------#
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


#--------------------------#
#     run STAN             #
#--------------------------#
stanFit <- stan(file = "orca_site_model.stan", 
                data = dataList, 
                pars = c("bage", "bpop", "bsex", "bind", "bsite", "bpopsite", "sigma", "permeth_pred"),
                warmup = 2000,
                iter = 12000, 
                chains = 3)


#-------------------------#
# Check MCMC Performance  #
#-------------------------#
print(stanFit)  
stan_trace(stanFit, pars = c("bage", "bsex", "bpop"), inc_warmup = TRUE)
stan_trace(stanFit, pars = "bind", inc_warmup = TRUE)
stan_trace(stanFit, pars = "bpopsite", inc_warmup = TRUE)
stan_trace(stanFit, pars = "sigma", inc_warmup = TRUE)


#-----------------------------------#
#  Plot Results With Stan Functions #
#-----------------------------------#
stan_plot(stanFit, par = c("bage", "bsex", "bpop"))
stan_plot(stanFit, par = "bind")
stan_plot(stanFit, par = "bsite")
stan_plot(stanFit, par = "sigma")
stan_plot(stanFit, par = "bpopsite")

s <- stan_plot(stanFit, par = c("bage", "bsex", "bpop"))
s + scale_y_discrete(labels=c("Age","Female","Male","Northern Resident","Southern Resident"))
s
str(s)

#############################################
#---     Posterior Predictive Check      ---#
#############################################

#-----------------------------------#
#      Extract Predicted Data       #         
#-----------------------------------#

#--- Extract the predictions ---#
mcmcChains = as.data.frame(stanFit)

write.csv(mcmcChains,'mcmcChains.csv')

chainLength = length(mcmcChains[, 1])

zypred = matrix(0, ncol = N, nrow = chainLength)
for (i in 1:N) {
  zypred[, i] = mcmcChains[, paste("permeth_pred[", i, "]", sep = "")]
}

#--- Mean expected value for record ---#
ypredMean = apply(zypred, 2, mean)

#--- Upper and lower expected 95% HDI for each visit ---#
ypredLow  = apply(zypred, 2, quantile, probs = 0.025)
ypredHigh = apply(zypred, 2, quantile, probs = 0.975)


#-----------------------------------#
#          Model Fit Plots          #         
#-----------------------------------#

#--- Plot mean predicted values ---#
par(mfrow = c(1, 1))
record = 1:500

subsample = sample(1:N, size = 500, replace=FALSE)



dotchart(ypredMean[subsample], xlim = c(-6, 6), xlab = "Standardized Percent Methylation", ylab = "Sample")
#--- Add HDIs ---#
segments(x0 = ypredLow[subsample], y0 = record, x1 = ypredHigh[subsample], y1 = record)
#--- Add observed values ---#
points(x = zpermeth[subsample], y = record, pch = 16, col = rgb(0, 0.8, 1, 0.6))


################################################################
#--------------------------------------------------------------#
#                PLOTS TO ASSESS EFFECTS                       #
#--------------------------------------------------------------#
################################################################


#-----------------------------------#
#  Plot Results With Stan Functions #
#          For Site Effects         #
#-----------------------------------#

actb<-c(which(startsWith(siteNames,"ACTB")))
bdnf<-c(which(startsWith(siteNames,"BDNF")))
crf<-c(which(startsWith(siteNames,"CRF")))
gapdh<-c(which(startsWith(siteNames,"GAPDH")))
ne<-c(which(startsWith(siteNames,"NE")))
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

gapdh_sites<-c()
for (i in gapdh){
  gapdh_sites<-append(gapdh_sites,paste("bsite[",i,"]", sep=""))
}

ne_sites<-c()
for (i in ne){
  ne_sites<-append(ne_sites,paste("bsite[",i,"]", sep=""))
}

np_sites<-c()
for (i in np){
  np_sites<-append(np_sites,paste("bsite[",i,"]", sep=""))
}

stan_plot(stanFit, par = c(actb_sites))
stan_plot(stanFit, par = c(gapdh_sites))
stan_plot(stanFit, par = c(bdnf_sites))
stan_plot(stanFit, par = c(crf_sites))
stan_plot(stanFit, par = c(ne_sites))
stan_plot(stanFit, par = c(np_sites))

#-----------------------------------#
#  Plot Results With Stan Functions #
#        For Pop*Site Effects       #
#-----------------------------------#


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

gapdh_ps<-c()
for (i in gapdh){
  gapdh_ps<-append(gapdh_ps,paste("bpopsite[1,",i,"]", sep=""))
  gapdh_ps<-append(gapdh_ps,paste("bpopsite[2,",i,"]", sep=""))
}

ne_ps<-c()
for (i in ne){
  ne_ps<-append(ne_ps,paste("bpopsite[1,",i,"]", sep=""))
  ne_ps<-append(ne_ps,paste("bpopsite[2,",i,"]", sep=""))
}

np_ps<-c()
for (i in np){
  np_ps<-append(np_ps,paste("bpopsite[1,",i,"]", sep=""))
  np_ps<-append(np_ps,paste("bpopsite[2,",i,"]", sep=""))
}

stan_plot(stanFit, par = c(actb_ps))
stan_plot(stanFit, par = c(gapdh_ps))
stan_plot(stanFit, par = c(bdnf_ps))
stan_plot(stanFit, par = c(crf_ps))
stan_plot(stanFit, par = c(ne_ps))
stan_plot(stanFit, par = c(np_ps))

plot(stanFit, par = c(actb_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(gapdh_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(bdnf_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(crf_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(ne_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")
plot(stanFit, par = c(np_ps), show_density = TRUE,
     ci_level = 0.95, fill_color = "skyblue")




#########################################################
### PLOT DIFFERENCES BETWEEN POPULATIONS AT EACH SITE ###
###                    SRKW - NRKW                    ###
#########################################################
mcmcChains = as.data.frame(stanFit)
chainLength = length(mcmcChains[, 1])

###---   ACTB   ---###

actb_ps_diff = c(rep(0,chainLength))
for (i in actb) {
  actb_ps_diff <- cbind(actb_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
}
actb_ps_diff<-as.data.frame(actb_ps_diff[,-1])
colnames(actb_ps_diff)<-c("-72","-79","-92","-94","-97","-108","-114","-117")
actb_diff<-pivot_longer(actb_ps_diff,1:length(actb_sites))

a_plot<-ggplot(actb_diff) +
       theme_bw(base_size = 22) +
       geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
       theme(legend.position = "none") +
       geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
       labs(y = "ACTB Sites") +
       theme(axis.title.x=element_blank())

actb_list <- c("-72", "-79", "-92", "-94", "-97", "-108", "-114", "-117")

a_plot<-ggplot(actb_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(actb_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "ACTB Sites") +
  theme(axis.title.x=element_blank())

###---   GAPDH   ---###
# 
# gapdh_ps_diff = c(rep(0,chainLength))
# for (i in gapdh) {
#   gapdh_ps_diff <- cbind(gapdh_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
# }
# gapdh_ps_diff<-as.data.frame(gapdh_ps_diff[,-1])
# colnames(gapdh_ps_diff)<-c(gapdh)
# gapdh_diff<-pivot_longer(gapdh_ps_diff,1:length(gapdh_sites))
# 
# g_plot<-ggplot(gapdh_diff) +
#       theme_bw(base_size = 22) +
#       geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
#       theme(legend.position = "none") +
#       geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
#       labs(y = "GAPDH Sites") +
#       theme(axis.title.x=element_blank(), axis.text.y=element_blank())
# 

###---   BDNF   ---###

bdnf_ps_diff = c(rep(0,chainLength))
for (i in bdnf) {
  bdnf_ps_diff <- cbind(bdnf_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
}
bdnf_ps_diff<-as.data.frame(bdnf_ps_diff[,-1])
colnames(bdnf_ps_diff)<-c("-70","-59","-51","-48","-45","-39","-15")
bdnf_diff<-pivot_longer(bdnf_ps_diff,1:length(bdnf_sites))


b_plot<-ggplot(bdnf_diff) +
      theme_bw(base_size = 22) +
      geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
      labs(y = "BDNF Sites") +
      theme(axis.title.x=element_blank())



# bdnf_list <- c("BDNF_26","BDNF_37","BDNF_45","BDNF_48","BDNF_51","BDNF_57","BDNF_81")
# b_plot<-ggplot(bdnf_diff) +
#  theme_bw(base_size = 22) +
#  geom_density_ridges(aes(x=value, y=factor(name, levels = c(bdnf_list)), alpha = 0.6), fill = "cyan") +
#  theme(legend.position = "none") +
#  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
#  labs(y = "BDNF Sites")+
#  theme(axis.title.x=element_blank())

###---   CRF   ---###

crf_ps_diff = c(rep(0,chainLength))
for (i in crf) {
  crf_ps_diff <- cbind(crf_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
}
crf_ps_diff<-as.data.frame(crf_ps_diff[,-1])
colnames(crf_ps_diff)<-c("-15", "-101", "-95", "-79", "-55", "-36", "-33")
crf_diff<-pivot_longer(crf_ps_diff,1:length(crf_sites))

c_plot<-ggplot(crf_diff) +
     theme_bw(base_size = 22) +
      geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
      labs(y = "CRF Sites")+
      theme(axis.title.x=element_blank())


crf_list <- c("-15", "-33", "-36", "-55", "-79", "-95", "-101")
c_plot<-ggplot(crf_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(crf_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "CRF Sites")+
  theme(axis.title.x=element_blank())

###---   NP   ---###

np_ps_diff = c(rep(0,chainLength))
for (i in np) {
  np_ps_diff <- cbind(np_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
}
np_ps_diff<-as.data.frame(np_ps_diff[,-1])
colnames(np_ps_diff)<-c("-34","-18","-15","-12","-9","-6","-1","-121","-100","-76","-74")
np_diff<-pivot_longer(np_ps_diff,1:length(np_sites))

np_plot<-ggplot(np_diff) +
      theme_bw(base_size = 22) +
      geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
      labs(y = "NR3C1 Sites") +
      theme(axis.title.x=element_blank())


np_list <- c("-1","-6","-9","-12","-15","-18","-34","-74","-76","-100","-121")
np_plot<-ggplot(np_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(np_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "NR3C1 Sites") +
  theme(axis.title.x=element_blank())

###---   NE   ---###
# 
# ne_ps_diff = c(rep(0,chainLength))
# for (i in ne) {
#   ne_ps_diff <- cbind(ne_ps_diff,(mcmcChains[,paste("bpopsite[2,",i,"]", sep="")]-mcmcChains[,paste("bpopsite[1,",i,"]", sep="")]))
# }
# ne_ps_diff<-as.data.frame(ne_ps_diff[,-1])
# colnames(ne_ps_diff)<-c(ne)
# ne_diff<-pivot_longer(ne_ps_diff,1:length(ne_sites))
# 
# ne_plot<-ggplot(ne_diff) +
#       theme_bw(base_size = 22) +
#       geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
#       theme(legend.position = "none") +
#       geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
#       labs(y = "NE Sites") +
#       theme(axis.title.x=element_blank(), axis.text.y=element_blank())

#################################################################
###                       6 PANEL PLOT                        ###
#################################################################


by_site_plots<-ggarrange(a_plot, g_plot, b_plot, c_plot, ne_plot, np_plot, nrow = 2, ncol = 3, labels = c("a","b","c","d","e","f"))

annotate_figure(by_site_plots, bottom = text_grob("Difference in Effect Size Between Populations", size = 24))

#################################################################
###                       4 PANEL PLOT                        ###
#################################################################


by_site_plots_4<-ggarrange(a_plot, b_plot, c_plot, np_plot, nrow = 2, ncol = 2, labels = c("A","B","C","D"))

annotate_figure(by_site_plots_4, bottom = text_grob("Difference in Effect Size Between Populations", size = 24))




#################################################################
### PLOT SUM OF DIFFERENCES BETWEEN POPULATIONS AT EACH LOCUS ###
###                      SRKW - NRKW                          ###
#################################################################


total_actb_diff<-data.frame(rowSums(actb_ps_diff), rep("ACTB",chainLength))
colnames(total_actb_diff)<-c("Difference","Locus")
total_gapdh_diff<-data.frame(rowSums(gapdh_ps_diff), rep("GAPDH",chainLength))
colnames(total_gapdh_diff)<-c("Difference","Locus")
total_bdnf_diff<-data.frame(rowSums(bdnf_ps_diff), rep("BDNF",chainLength))
colnames(total_bdnf_diff)<-c("Difference","Locus")
total_crf_diff<-data.frame(rowSums(crf_ps_diff), rep("CRF",chainLength))
colnames(total_crf_diff)<-c("Difference","Locus")
total_ne_diff<-data.frame(rowSums(ne_ps_diff), rep("NE",chainLength))
colnames(total_ne_diff)<-c("Difference","Locus")
total_np_diff<-data.frame(rowSums(np_ps_diff), rep("NR3C1",chainLength))
colnames(total_np_diff)<-c("Difference","Locus")

total_all_diff<-bind_rows(total_actb_diff, total_gapdh_diff, total_bdnf_diff, total_crf_diff, total_ne_diff, total_np_diff)


ggplot(total_all_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=Difference, y=factor(Locus, levels = c("NR3C1","NE","CRF","BDNF","GAPDH","ACTB")), alpha = 0.6), fill="salmon") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(x = "Difference in Effect Between Populations Across All Sites") +
  theme(axis.title.y=element_blank()) + 
  xlim(-6,6)

###  FOUR LOCI FOR FIGURES  ###

total_four_diff<-bind_rows(total_actb_diff, total_bdnf_diff, total_crf_diff, total_np_diff)


ggplot(total_four_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=Difference, y=factor(Locus, levels = c("NR3C1","CRF","BDNF","ACTB")), alpha = 0.6), fill="salmon") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(x = "Difference in Effect Between Populations Across All Sites") +
  theme(axis.title.y=element_blank()) + 
  xlim(-6,6)
  
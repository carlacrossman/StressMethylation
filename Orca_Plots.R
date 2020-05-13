#########################################################
### PLOT DIFFERENCES BETWEEN POPULATIONS AT EACH SITE ###
###                    SRKW - NRKW                    ###
#########################################################

### --- Load Libraries --- ###
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("C:/Users/c_cro/Documents/PhD/Epigenetics/Methylation/results")
options(mc.cores = parallel::detectCores())
mcmcChains = read.csv("mcmcChains.csv")
chainLength = length(mcmcChains[, 1])

#####################################################
###   Import the original data file to pull out   ###
###      demographic or other data as needed      ###
#####################################################

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
siteNames = levels(orca$site)


### --- Identify the list of sites for each locus --- ###
actb<-c(which(startsWith(siteNames,"ACTB")))
bdnf<-c(which(startsWith(siteNames,"BDNF")))
crf<-c(which(startsWith(siteNames,"CRF")))
np<-c(which(startsWith(siteNames,"NP")))


#################################################################
###              Fig 2: NO SIG & ACTB PLOT                    ###
#################################################################

### --- SETTING UP DATA STRUCTURE TO BE READ INTO GGPLOT RIDGES OBJECT --- ###

age = mcmcChains[, "bage"]
NRKW = mcmcChains[, "bpop.1."]
SRKW = mcmcChains[, "bpop.2."]
female = mcmcChains[, "bsex.1."]
male = mcmcChains[, "bsex.2."]


age.frame = data.frame(rep("Age", length(age)), age)
colnames(age.frame) = c("variable","posteriors")
NRKW.frame = data.frame(rep("NRKW", length(NRKW)), NRKW)
colnames(NRKW.frame) = c("variable","posteriors")
SRKW.frame = data.frame(rep("SRKW", length(SRKW)), SRKW)
colnames(SRKW.frame) = c("variable","posteriors")
female.frame = data.frame(rep("Female", length(female)), female)
colnames(female.frame) = c("variable","posteriors")
male.frame = data.frame(rep("Male", length(male)), male)
colnames(male.frame) = c("variable","posteriors")

demographics = rbind(age.frame, NRKW.frame, SRKW.frame, female.frame, male.frame)

variable_list = c("Male", "Female", "SRKW", "NRKW", "Age")


### --- Setting up code for the posteriors of age, population and sex --- ###
demo<-ggplot(demographics) +
        theme_bw(base_size = 18) +
        geom_density_ridges(aes(x=posteriors, y=factor(variable, levels = c(variable_list)), alpha = 0.6), fill = "cyan") +
        theme(legend.position = "none") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
        labs(x = "Posterior Probability of Effects") +
        theme(axis.title.y=element_blank())


########################### PART B OF THIS FIGURE FOR ACTB ###########################

### --- Pull out ACTB Sites --- ###
actb_sites<-c()
for (i in actb){
  actb_sites<-append(actb_sites,paste("bsite[",i,"]", sep=""))
}

### --- Calculate the difference in posteriors between the two populations --- ###
actb_ps_diff = c(rep(0,chainLength))
for (i in actb) {
  actb_ps_diff <- cbind(actb_ps_diff,(mcmcChains[,paste("bpopsite.2.",i,".", sep="")]-mcmcChains[,paste("bpopsite.1.",i,".", sep="")]))
}
actb_ps_diff<-as.data.frame(actb_ps_diff[,-1])
colnames(actb_ps_diff)<-c("-72","-79","-92","-94","-97","-108","-114","-117")
actb_diff<-pivot_longer(actb_ps_diff,1:length(actb_sites))


### --- Generate the ACTB difference in posteriors plot --- ###
a_plot<-ggplot(actb_diff) +
  theme_bw(base_size = 22) +
  geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "ACTB Sites") +
  theme(axis.title.x=element_blank())

### --- Reorder the sites for the figure --- ###
actb_list <- c("-72", "-79", "-92", "-94", "-97", "-108", "-114", "-117")

### --- Plot posteriors for ACTB --- ###
a_plot<-ggplot(actb_diff) +
  theme_bw(base_size = 18) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(actb_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "ACTB Sites", x = "Difference in Posterior Probability of Effects")


#####################  ASSEMBLE FIGURE ###########################

### --- Combine demogrphic and ACTB into a 2 panel figure --- ###
demo_actb<-ggarrange(demo, a_plot, nrow = 1, ncol = 2, labels = c("A","B"))
demo_actb


#################################################################
###               Figure 3 - 3 PANEL PLOT                     ###
###     Difference at CpG sites in 3 stress related loci      ###
#################################################################

###########################   BDNF   ###########################
### --- Identify BDNF sites --- ###
bdnf_sites<-c()
for (i in bdnf){
  bdnf_sites<-append(bdnf_sites,paste("bsite[",i,"]", sep=""))
}

### --- Calculate the in posteriors between populations --- ###
bdnf_ps_diff = c(rep(0,chainLength))
for (i in bdnf) {
  bdnf_ps_diff <- cbind(bdnf_ps_diff,(mcmcChains[,paste("bpopsite.2.",i,".", sep="")]-mcmcChains[,paste("bpopsite.1.",i,".", sep="")]))
}
bdnf_ps_diff<-as.data.frame(bdnf_ps_diff[,-1])
colnames(bdnf_ps_diff)<-c("-70","-59","-51","-48","-45","-39","-15")
bdnf_diff<-pivot_longer(bdnf_ps_diff,1:length(bdnf_sites))

### --- Generate BDNF posterior plot --- ###
b_plot<-ggplot(bdnf_diff) +
  theme_bw(base_size = 18) +
  geom_density_ridges(aes(x=value, y=name, alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "BDNF Sites") +
  theme(axis.title.x=element_blank())



###########################   CRF   ###########################
### --- Identify CRF sites --- ###
crf_sites<-c()
for (i in crf){
  crf_sites<-append(crf_sites,paste("bsite[",i,"]", sep=""))
}

### --- Calculate the in posteriors between populations --- ###
crf_ps_diff = c(rep(0,chainLength))
for (i in crf) {
  crf_ps_diff <- cbind(crf_ps_diff,(mcmcChains[,paste("bpopsite.2.",i,".", sep="")]-mcmcChains[,paste("bpopsite.1.",i,".", sep="")]))
}
crf_ps_diff<-as.data.frame(crf_ps_diff[,-1])
colnames(crf_ps_diff)<-c("-15", "-101", "-95", "-79", "-55", "-36", "-33")
crf_diff<-pivot_longer(crf_ps_diff,1:length(crf_sites))

### --- Generate CRF posterior plot --- ###
crf_list <- c("-15", "-33", "-36", "-55", "-79", "-95", "-101")
c_plot<-ggplot(crf_diff) +
  theme_bw(base_size = 18) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(crf_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "CRF Sites")+
  theme(axis.title.x=element_blank())

###########################   NP   ###########################
### --- Identify NR3C1 sites --- ###
np_sites<-c()
for (i in np){
  np_sites<-append(np_sites,paste("bsite[",i,"]", sep=""))
}

### --- Calculate the in posteriors between populations --- ###
np_ps_diff = c(rep(0,chainLength))
for (i in np) {
  np_ps_diff <- cbind(np_ps_diff,(mcmcChains[,paste("bpopsite.2.",i,".", sep="")]-mcmcChains[,paste("bpopsite.1.",i,".", sep="")]))
}
np_ps_diff<-as.data.frame(np_ps_diff[,-1])
colnames(np_ps_diff)<-c("-34","-18","-15","-12","-9","-6","-1","-121","-100","-76","-74")
np_diff<-pivot_longer(np_ps_diff,1:length(np_sites))

### --- Generate NR3C1 posterior plot --- ###
np_list <- c("-1","-6","-9","-12","-15","-18","-34","-74","-76","-100","-121")
np_plot<-ggplot(np_diff) +
  theme_bw(base_size = 18) +
  geom_density_ridges(aes(x=value, y=factor(name, levels = c(np_list)), alpha = 0.6), fill = "cyan") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", size = 2) +
  labs(y = "NR3C1 Sites") +
  theme(axis.title.x=element_blank())


#####################  ASSEMBLE FIGURE ###########################

### --- Compile plots into a three panel figure --- ###
by_site_plots_3<-ggarrange(b_plot, np_plot, c_plot, nrow = 1, ncol = 3, labels = c("A","B","C"))
annotate_figure(by_site_plots_3, bottom = text_grob("Difference in Posterior Probability of Effects Between Populations", size = 20))

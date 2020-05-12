############################################################################
# 																                                 	       #
#   IDENTIFYING C/T RATIOS AT CG SITE AND PREPARING DATA FRAME FOR MODEL   #
#		NOTE SEQUENCING RUN CONTAINED MORE GENES THAN WERE INCLUDED IN THIS    #
#   STUDY. THESE GENES AND SITES ARE EXCLUDED WITH THE FIRST LINES OF      #
#   CODE IN THE BAYESIAN MODEL FILE : Methylation_model_4_loci.R           #
# 																                                 	       #
############################################################################
### --- LOAD LIBRARIES --- ###

library(BiocManager)   # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(Biostrings)    # installed with above as BiocManager::install("Biostrings")
library(Rsamtools)     # installed with above as BiocManager::install("Rsamtools")
library(ShortRead)     # installed with above as BiocManager::install("ShortRead")
library(GenomicRanges) 
library(ape)
library(ggplot2)
library(stringr)


### --- SET PARAMETERS FOR REFERENCE SEQUENCES --- ###
#------    This helps set positions from the   ------#
#------    reference seq that can be refered   ------# 
#------         to in the other sequences      ------# 

params<-data.frame(c("BDNF","CRF","NP","ACTB","GAPDH","NE"),c("ref_BDNF:1-276","ref_CRF:1-351","ref_NP:1-232","ref_ACTB:1-163","ref_GAPDH:1-183","ref_NE:1-335"))   # Create a data frame with loci and ref_seq name:length for indexing
colnames(params)<-c("locus","param")                             # Set column names for parameter data frame


### --- COMPILE LIST OF CG SITES --- ###

BDNF_CG<-c(26,37,45,48,51,57,81,106,145,154,167,172,192,208,213)                                # Provide a vector of CG sites by poistion in the BDNF reference sequence
CRF_CG<-c(27,33,49,73,92,95,113,135,148,236,242,257,260,265,270,286,305,313,319,321)            # Provide a vector of CG sites by poistion in the CRF reference sequence
NP_CG<-c(45,66,90,92,132,148,151,154,157,160,165,176,196,201)                                   # Provide a vector of CG sites by poistion in the GAPDH reference sequence
GAPDH_CG<-c(35,55,58,89,117,126,131)                                                            # Provide a vector of CG sites by poistion in the NE reference sequence
NE_CG<-c(1,7,25,29,44,53,58,64,70,72,83,85,96,101,103,108,126,136,144,150,160,168,172,179,189,196,214,227,250,261,268,287,294,297,300)            # Provide a vector of CG sites by poistion in the NP reference sequence
ACTB_CG<-c(3,24,29,36,49,51,54,65,71,74,162)                                                    # Provide a vector of CG sites by poistion in the ACTB reference sequence

CG_site_list<-list(BDNF=BDNF_CG,CRF=CRF_CG,NP=NP_CG,ACTB=ACTB_CG,GAPDH=GAPDH_CG,NE=NE_CG)       # Create a list of CG sites for all loci. Each list element is a vector of CG sites for a locus. 


######################################################################
#       ITERATE TRHOUGH EACH SAMPLE TO GENERATE A DATA FRAME         #
#           WITH THE NUMBER OF C & T COUNTS AT CG SITES              #
######################################################################


CG_sites=data.frame(matrix(ncol = 5, nrow = 0))                                      # Create an empty data frame results will be written to

for (locus in c("BDNF","CRF","NP","GAPDH","NE","ACTB")){                             # Iterate through all loci

  path = paste("C:/Users/c_cro/Desktop/basespace/sorted/",locus,sep="")              # Set Directory
  file.names<-dir(path,pattern=".bam$")                                              # Identify unique files
  
  for (sample in file.names) {                                                       # Open for loop to iterate through all samples
    locus_bam<-BamFile(paste(path,sample,sep="/"))                                   # Load Bam Files
    stringbam_locus<-stackStringsFromBam(locus_bam, param=as.character(params[params$locus==locus,2]))      # Stack bam file to create a list of sequence reads
    seqs<-DNAStringSet(stringbam_locus)                                              # Specify the list of sequences as DNA and write to a file "seqs"
    
    for (site in CG_site_list[[locus]]) {                                            # Open for loop to iterate through known CpG sites specified above
      Base_Count<-as.data.frame(table(subseq(seqs, start=site, end=site, width=NA))) # Pull out the base counts for the desired site and write to a data frame
      C_Bases<-Base_Count[Base_Count[,1]=="C",2]                                     # Write the number of C Bases at this site to a temporary variable
      T_Bases<-Base_Count[Base_Count[,1]=="T",2]                                     # Write the number of T Bases at this site to a temporary variable
      if (length(C_Bases) == 0){C_Bases<-0}                                          # If there were no C Bases, write 0 to vector
      if (length(T_Bases) == 0){T_Bases<-0}                                          # If there were no T Bases, write 0 to vector
      g.site<-paste(locus,site,sep="_")                                              # Create a site name that merges current locus and current site
      newrow<-data.frame(sample, locus, g.site, C_Bases, T_Bases)                    # Write a new data frame with sample no, locus, CG site, # Cs and # Ts
      CG_sites<-rbind(CG_sites, newrow)                                              # Append this new data to an existing data frame
    }                                                                                # Close for loop for iterating through CpG sites
  }                                                                                  # Close for loop for iterating through samples
}  

######################################################################
#       PROVIDE THE CG_sites DATA FRAME WITH MORE INFORMATION        #
######################################################################

### --- MANIPULATE THE COLLECTED DATA TO PROVIDE MEANINGFUL C_RATIOS --- ###

colnames(CG_sites)=c("well","locus","site","C","T")                                # Rename columns for data frame with C and T data
CG_sites$"total"   <- CG_sites$C + CG_sites$T                                      # Add column for total number of C or T calls for the CG sites 
CG_sites$"C_Ratio" <- CG_sites$C / CG_sites$total                                  # Add column that calculates the C/T ratio which is a metric for percent methylation


### --- MANIPULATE THE WELL NAMES TO PROVIDE MEANINGFUL INDEXES --- ###

CG_sites$index<-c()                                                                # Create an empty column for index data will be stored
for (i in 1:nrow(CG_sites)){                                                       # Iterate through each row of the data frame
  CG_sites$index[i]<-str_split(as.character(CG_sites$well[i]),"_")[[1]][2]         # To each row of the index, write the well number by splitting the character string provided in the well column and retaining the plate location (well number) as the index
}                                                                                  # Close for loop iterating through each row of the data frame


### --- IMPORT PLATE LOCATION SPREADSHEET AND MERGE TO CG_sites --- ###

plate<-read.csv("KWdatabase_plate.csv")                                            # Import plate information. All loci at each locus should have a sample number even if NA
colnames(plate)<-c("Sample.No","locus","index")                                    # Rename columns for plate data frame 
m<-merge(CG_sites,plate,by=c("locus","index"),all=TRUE)                            # Merge Sample No. into CG sites data frame from plate by locus and index


### --- IMPORT SAMPLE SPREADSHEET AND MERGE TO CG_sites --- ###

samples<-read.csv("KWdatabase_sample.csv")                                         # Import sample information.
samp<-samples[,c(1,2,4)]                                                           # Extract relevant columns from sample file to just use Sample.No, ID and age
colnames(samp)<-c("Sample.No","ID","age")                                          # Rename columns for sample data frame
m2<-merge(m,samp, by="Sample.No")                                                  # Merge age and ID into CG_site data frame by Sample No.


### --- IMPORT INDIVIDUAL SPREADSHEET AND MERGE TO CG_sites --- ###

ind<-read.csv("KWdatabase_individual.csv")                                         # Import individual information.
ind2<-ind[,1:3]                                                                    # Extract relevant columns from sample file to just use ID, sex and population
CG_Sites_Info<-merge(m2,ind2, by="ID")                                             # Merge sex and population into CG_site data frame by ID


### --- INDENTIFY READ DEPTH THRESHOLDS FOR INCLUSION --- ###
sub_BDNF<-CG_Sites_Info[CG_Sites_Info$locus=="BDNF",]
hist(sub_BDNF$total)
# 5000 for BDNF

sub_CRF<-CG_Sites_Info[CG_Sites_Info$locus=="CRF",]
hist(sub_CRF$total)
# 100 for CRF

sub_GAPDH<-CG_Sites_Info[CG_Sites_Info$locus=="GAPDH",]
hist(sub_GAPDH$total)
# 20000 for GAPDH

sub_ACTB<-CG_Sites_Info[CG_Sites_Info$locus=="ACTB",]
hist(sub_ACTB$total)
# 10000 for ACTB

sub_NE<-CG_Sites_Info[CG_Sites_Info$locus=="NE",]
hist(sub_NE$total)
# 100 for NE

sub_NP<-CG_Sites_Info[CG_Sites_Info$locus=="NP",]
hist(sub_NP$total)
# 500 for NP

## Remove rows that don't meet minimum read depth thresholds
too_short<-c(which(CG_Sites_Info$total<500 & CG_Sites_Info$locus=="NP"),
which(CG_Sites_Info$total<100 & CG_Sites_Info$locus=="NE"),
which(CG_Sites_Info$total<10000 & CG_Sites_Info$locus=="ACTB"),
which(CG_Sites_Info$total<20000 & CG_Sites_Info$locus=="GAPDH"),
which(CG_Sites_Info$total<100 & CG_Sites_Info$locus=="CRF"),
which(CG_Sites_Info$total<5000 & CG_Sites_Info$locus=="BDNF"))

CG_Sites_Info_Trimmed<-CG_Sites_Info[-too_short,]

### --- WRITE CG_SITES TO OUTOUT FILE --- ###

write.csv(CG_Sites_Info_Trimmed, file="orcaData.csv")



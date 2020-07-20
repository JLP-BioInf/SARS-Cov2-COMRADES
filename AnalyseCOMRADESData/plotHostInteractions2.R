################################################################
# This code is specific to the project:                        #
# The short- and long-range RNA-RNA Interactome of SARS-CoV-2  #
#  It uses output from plotHostInteractions1.R to or other     #
#   RNAs of interest to plot where it binds along the genome   #
################################################################
#author -  Jonathan Price

source(COMRADESfunctions.R)





##############################
# Make libraries 
##############################

library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(igraph)
library(heatmap3)
library(MASS)
library(mixtools)



##############################
# Set Variables
##############################

hybDirSars = "../data/sars/"
IDSars  = "NC045512_hsapiens_SARSCOV-1_SARSCOV"


hybDirMers = "../data/mers/"
IDMers  = "NC038294_hsapians_MERS-1_MERS"






##############################
# Load in the Files 
##############################
hybListMers = readHybFiles(hybDirMers)
hybListSars = readHybFiles(hybDirSars)





##############################
# snRNA
##############################


# Choose the features of interest
featureList = list()
featureList[[IDMers]] = read.table(text = "ENSCSAT00000021419_ENSCSAG00000021329_U4_snRNA
ENSCSAT00000027369_ENSCSAG00000027279_U4_snRNA
ENSCSAT00000026089_ENSCSAG00000025999_U2_snRNA
ENSCSAT00000020532_ENSCSAG00000020442_U2_snRNA
ENSCSAT00000026292_ENSCSAG00000026202_U1_snRNA
ENSCSAT00000022634_ENSCSAG00000022544_U1_snRNA
ENSCSAT00000025579_ENSCSAG00000025489_U1_snRNA
ENSCSAT00000028010_ENSCSAG00000027920_U1_snRNA
ENSCSAT00000023365_ENSCSAG00000023275_U1_snRNA
ENSCSAT00000028032_ENSCSAG00000027942_U1_snRNA
ENSCSAT00000023715_ENSCSAG00000023625_U1_snRNA
ENSCSAT00000025680_ENSCSAG00000025590_U1_snRNA
ENSCSAT00000026933_ENSCSAG00000026843_U1_snRNA
ENSCSAT00000022422_ENSCSAG00000022332_U1_snRNA
ENSCSAT00000027400_ENSCSAG00000027310_U1_snRNA
ENSCSAT00000020657_ENSCSAG00000020567_U1_snRNA
ENSCSAT00000025955_ENSCSAG00000025865_U1_snRNA
ENSCSAT00000020456_ENSCSAG00000020366_U1_snRNA", header = F)$V1


#genome length
lengthList = list()
lengthList[[IDMers]] = c(0,30000)



for(TE in IDMers){
  
  
  
  dataF = data.frame()
  
  for(feature in featureList[[TE]]){
    
    

    for(hyb in 1:length(hybListMers)){
      #get the vectors of interaction sites for the different samples
      sampleVect = getVectorofInteractions(TE, feature, hybListMers[[hyb]])
      sampleVect$id = "sample"
      if(hyb %% 2 != 0 ){
        sampleVect$id = "control"
      }
      sampleVect$feature = feature
      if(!(is.data.frame(sampleVect))){next}
      dataF = rbind(dataF, sampleVect)
    }
  }
  
  
  dataF = aggregate(dataF$Freq, by= list(dataF$vect,dataF$id, dataF$feature), FUN = sum)
  colnames(dataF) = c("vect","id","feature","Freq")
  dataF$featureReduced = str_split_fixed(dataF$feature,"_",4)[,3]

  
  pdf(paste("../figures/specificInteractionGraphs/",TE,"_7sk_Combined_all.pdf", sep = ""), height = 2, width = 9 )
  max = round(max(as.numeric(as.character(dataF$vect))),-(nchar(max(as.numeric(as.character(dataF$vect))))-1))
  plot(ggplot()+
         geom_vline(xintercept = seq(0,max,1000), colour = "grey", alpha = 0.6)+
         geom_vline(xintercept = seq(0,max,500), colour = "grey", alpha = 0.4, linetype = "dashed")+
         geom_vline(xintercept = seq(0,max,250), colour = "grey", alpha = 0.3, linetype = "dotted")+
         geom_bar(data = dataF,aes(x = as.numeric(as.character(vect)), y = as.numeric(as.character(Freq))), stat = "identity")+
         ggtitle(paste(TE,"and",feature))+
         facet_grid(id~., scales = "free_y")+
         theme_classic())+
    theme(legend.position = "none")
  dev.off()
  
  
  pdf(paste("../figures/specificInteractionGraphs/",TE,"_7sk_split.pdf", sep = ""), height = 3, width = 9 )
  plot(ggplot()+
         geom_vline(xintercept = seq(0,max,1000), colour = "grey", alpha = 0.6)+
         geom_vline(xintercept = seq(0,max,500), colour = "grey", alpha = 0.4, linetype = "dashed")+
         geom_vline(xintercept = seq(0,max,100), colour = "grey", alpha = 0.3, linetype = "dotted")+
         geom_bar(data = dataF,aes(x = as.numeric(as.character(vect)), y = as.numeric(as.character(Freq)), fill = id), stat = "identity")+
         ggtitle(paste(TE,"and",feature))+
         facet_grid(featureReduced+id~.)+
         theme_classic())+
    theme(legend.position = "none")
  dev.off()
  
  
  
  
}




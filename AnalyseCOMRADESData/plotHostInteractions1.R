################################################################
# This code is specific to the project:                        #
# The short- and long-range RNA-RNA Interactome of SARS-CoV-2  #
#  It processes the chimeric reads to find important host      #
#   interactions                                               #
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



hybDirSarsmRNA = "../data/sarsmRNA/"


##############################
# Load in the Files 
##############################
hybListMers = readHybFiles(hybDirMers)
hybListSars = readHybFiles(hybDirSars)


hybListSarsmRNA = readHybFiles(hybDirSarsmRNA)
hybList = hybListSarsmRNA

#find ostr abundat mrna
head(table(hybList[[2]]$V4)[order(table(hybList[[2]]$V4), decreasing = T)])

##############################
# From each file parse so that there are 2 columns each with the TE in column 1 and the interactor in column 2
##############################
#Just grab the columns relating to the feature that is interacting
TE = "NC045512_hsapiens_SARSCOV-1_N"
alteredHybList = list()
for(hyb in 1:length(hybList)){
  
  controlHyb = hybList[[hyb]]
  #get the right columns
  controlHybO  = as.data.frame(cbind(as.character(controlHyb$V4), as.character(controlHyb$V10)))
  #zika as ref
  TEs = c(TE)
  for(TE in TE){
    # alteredHybList[[TE]] = list()
    #Only get those lines with the zika in
    controlHyb = controlHybO[controlHybO$V1 == TE | controlHybO$V2 == TE,]
    
    #remove the factor levels 
    controlHyb$V1 = as.character(controlHyb$V1)
    controlHyb$V2 = as.character(controlHyb$V2)
    
    
    #Swap the columns around that do not have zika in column 2 
    controlHybtmp1 = controlHyb[controlHyb$V1 == TE,]
    controlHybtmp2 = controlHyb[!(controlHyb$V1 == TE),]
    controlHybtmp2 = rev(controlHybtmp2)
    colnames(controlHybtmp2) = c("V1","V2")
    controlHyb = as.data.frame(rbind(controlHybtmp1, controlHybtmp2))
    
    # check to see if the unique stuff has reduced 
    print(unique(sort(controlHyb$V1)))
    print(length(unique(sort(controlHyb$V2))))
    
    
    #add to list
    alteredHybList[[TE]][[hyb]] = controlHyb
  }
}




##############################
# Now Combine the two  
##############################

hybDir = hybDirSarsmRNA
for(TE in TE){
  aggList = list()
  totalNames = c()
  for(hyb in 1:length(alteredHybList[[TE]])){
    #aggList[[TE]] = list()
    sampleHyb = alteredHybList[[TE]][[hyb]]
    sData = sampleHyb[sampleHyb$V1 == TE,]
    freqSample2 = aggregate(sData$V1, by = list(sData$V2), FUN = length)
    freqSample = freqSample2$x
    names(freqSample) = freqSample2$Group.1
    aggList[[TE]][[hyb]] = freqSample
    # Get the total features that exist in the dataset 
    totalNames = unique(sort(c(totalNames, names(freqSample))))
  }
  
  tmpMat = matrix(0, nrow = length(totalNames), ncol = length(hybList))
  row.names(tmpMat) = totalNames
  colnames(tmpMat) = list.files(hybDir, pattern = "hybrids.hyb$")
  
  
  
  for(i in 1:nrow(tmpMat)){
    for(j in 1:ncol(tmpMat)){
      tmpMat[i,j] = aggList[[TE]][[j]][row.names(tmpMat)[i]]
    }
    print(i)
  }
  
  tmpMat[is.na(tmpMat)] <- 0.00001
  write.table(tmpMat, file = paste("../hybStats2/TEs_",TE,"mRNA_featureStatsNotMerged.txt",sep = ""), quote = F, sep = "\t")
}















library(reshape2)

list.files("../hybStats2/",pattern = "*mRNA*.*txt$",full.names = T)



statsList = list()

for(i in 1){
  statsList[[i]] = read.table( list.files("../hybStats2/",pattern = "*mRNA*.*txt$",full.names = T)[i], header = T, row.names = 1)
  statsList[[i]]$ID = sapply(row.names(statsList[[i]]), function(x) strsplit(x, "_")[[1]][4], USE.NAMES=FALSE)

}
statsList 

# Now get the stats from aggregate 

aggList = list()
aggList2 = list()

for(i in 1){
  
  aggList[[i]] =   aggregate(.~ID , statsList[[i]], sum)

  
  factor1 = sum(aggList[[i]]$SARS2m_1S.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb) /sum(aggList[[i]]$SARS2m_1C.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb)
  factor2 = sum(aggList[[i]]$SARS2m_2S.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb) /sum(aggList[[i]]$SARS2m_2C.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb)
  
  
  aggList[[i]]$M1Norm = aggList[[i]]$SARS2m_1S.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb/ (aggList[[i]]$SARS2m_1C.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb * factor1)
  aggList[[i]]$M2Norm = aggList[[i]]$SARS2m_2S.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb/ (aggList[[i]]$SARS2m_1C.assembled.tstk_comp_ChlSab1_and_SARS_plus_mRNA_hybrids.hyb* factor2)
  data = melt(aggList[[i]],id.vars = c("ID") )
  print(data)

  
  pdf(paste("../hybStats2/stats_",list.files("../hybStats2/",pattern = "*mRNA*.*txt$")[i], ".pdf", sep = ""), height = 4, width = 7)
  
  plot(ggplot() + 
         geom_boxplot(data = data[data$variable %in% c("M1Norm","M3Norm", "M2Norm"),], aes(x = reorder(ID, value, FUN = median), y = log2(value))) +
         geom_point(data = data[data$variable %in% c("M1Norm","M3Norm", "M2Norm"),], aes(x = reorder(ID, value, FUN = median), y = log2(value), fill = ID)) +
         geom_hline(yintercept = 0, colour = "darkred")+
         theme_classic() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
  
}












aggList = list()
aggList2 = list()
for(i in 1:length(statsList)){
  
  aggList[[i]] =   aggregate(.~ID , statsList[[i]], sum)
  
  factor1 = sum(aggList[[i]][,10])/sum(aggList[[i]][,9])
  
  
  aggList[[i]]$P1Norm = aggList[[i]]$newPolysomes / (aggList[[i]]$newPolsomes.cont * factor1)
  
  data = melt(aggList[[i]],id.vars = c("ID") )
  print(data)
  data = data[data$variable %in% c("P1Norm") & data$ID != "Ig" & data$ID != "microRNA" & data$ID != "snoRNA" &
                data$ID != "snRNA" & data$ID != "tRNA" & data$ID != "miscRNA" & data$ID != "microRNA",]
  
  pdf(paste("stats_",list.files(pattern = ".*txt$")[i], ".pdf", sep = ""), height = 3, width = 4)
  
  plot(ggplot() + 
         geom_boxplot(data = data, aes(x = reorder(ID, value, FUN = median), y = log2(value))) +
         geom_point(data = data, aes(x = reorder(ID, value, FUN = median), y = log2(value), fill = ID)) +
         geom_hline(yintercept = 0, colour = "darkred")+
         theme_classic() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  dev.off()
  
  plot(    ggplot() + 
             geom_boxplot(data = data, aes(x = reorder(ID, value, FUN = median), y = log2(value))) +
             geom_point(data = data, aes(x = reorder(ID, value, FUN = median), y = log2(value), fill = ID)) +
             geom_hline(yintercept = 0, colour = "darkred")+
             theme_classic() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
}








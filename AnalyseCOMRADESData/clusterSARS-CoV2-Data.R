################################################################
# This code is specific to the project:                        #
# The short- and long-range RNA-RNA Interactome of SARS-CoV-2  #
#  It clusters the SARS Hyb data to create chimeric groups     #
################################################################
#author -  Jonathan Price



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
library(RColorBrewer)
library(foreach)
library(doParallel)
source(COMRADESfunctions.R)

##############################
# Set Variables
##############################

hybDirSars = "../data/sars/"
IDSars  = "NC045512_hsapiens_SARSCOV-1_SARSCOV"


##############################
# Load in the Files 
##############################

hybListSars = readHybFiles(hybDirSars)


##############################
# Make a matrix for the contact maps
# and plot (inlcuding the line traces )
##############################



hybListSarsSwap = swapHybs(hybListSars, IDSars)
#matrixListSars = getMatrices(hybListSarsSwap, IDSars)






##############################
# clustering 
##############################

##########

longDistHyb = subsetHybList(hybListSarsSwap,10,30000,length = 800)



##############################
# From these long range interactions 
#find groups of reads. 
##############################


##############################
# Make the Granges
##############################
chimeraList = hybToGRanges(longDistHyb,IDSars)




##############################
# Redduce the GRanges to make them smaller
##############################



chimeraListSampled =list()

  for(i in 1:length(chimeraList)){
    max =  length(chimeraList[[i]][["left"]])
    seq = c(1,max)
    if(max > 40000){seq = seq(1,max,by = 20000)
    seq = c(seq,max)}
    
    print(seq)
    chimeraListSampled[[i]] = list()
    for(j in c("left","right","gap")){
      chimeraListSampled[[i]][[j]] = list()
      for(k in 1:(length(seq)-1)){
        sample = seq[k]:seq[k+1]
        print(seq[k])
        print(seq[k+1])
      chimeraListSampled[[i]][[j]][[k]] = chimeraList[[i]][[j]][sample]
    }
  }
  
}



registerDoParallel(5) 

plottingList = list()
for(i in 1:4){

  plottingList[[i]] = list()
foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
  adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],"nucleotide", 15)
  net = graph_from_adjacency_matrix(adjacancyMat, mode = "undirected", weighted = T)
  clustering = cluster_walktrap(net,steps = 2)
  highest_clusters = names(table(membership(clustering))[table(membership(clustering)) > 10])
  plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters, 
                                              chimeraListSampled[[i]][["left"]][[k]], 
                                              chimeraListSampled[[i]][["right"]][[k]])
  
  print(k)
}

}
#adjacancyMat = getAdjacancyMat(chimeraListSampled[[2]][["gap"]],"perc", 0.7)

longRange = plottingList
#save(longRange, file = "./longRangeclustering.RData")

##############################
# Now do the clustering for short range
##############################


longDistHyb = subsetHybList(hybListSarsSwap,1,9,length = 800)




##############################
# Make the Granges
##############################
chimeraList = hybToGRanges(longDistHyb,IDSars)




##############################
# Redduce the GRanges to make them smaller
##############################



chimeraListSampled =list()

for(i in 1:length(chimeraList)){
  max =  length(chimeraList[[i]][["left"]])
  seq = c(1,max)
  if(max > 40000){seq = seq(1,max,by = 20000)
  seq = c(seq,max)}
  
  print(seq)
  chimeraListSampled[[i]] = list()
  for(j in c("left","right","gap")){
    chimeraListSampled[[i]][[j]] = list()
    for(k in 1:(length(seq)-1)){
      sample = seq[k]:seq[k+1]
      print(seq[k])
      print(seq[k+1])
      chimeraListSampled[[i]][[j]][[k]] = chimeraList[[i]][[j]][sample]
    }
  }
  
}



registerDoParallel(5) 

plottingList = list()
for(i in 1:4){
  
  plottingList[[i]] = list()
  foreach (k=1:length( chimeraListSampled[[i]][["gap"]])) %do% {
    adjacancyMat = getAdjacancyMat(chimeraListSampled[[i]][["gap"]][[k]],"nucleotide", 15)
    net = graph_from_adjacency_matrix(adjacancyMat, mode = "undirected", weighted = T)
    clustering = cluster_walktrap(net,steps = 2)
    highest_clusters = names(table(membership(clustering))[table(membership(clustering)) > 10])
    plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters, 
                                                chimeraListSampled[[i]][["left"]][[k]], 
                                                chimeraListSampled[[i]][["right"]][[k]])
    
    print(k)
  }
  
}
#adjacancyMat = getAdjacancyMat(chimeraListSampled[[2]][["gap"]],"perc", 0.7)

shortRange = plottingList
#save(shortRange, file = "./shortRangeclustering.RData")





##############################
# Add the short and the long range together
##############################




combinedPlotting = shortRange
for(i in 1:length(shortRange)){
  plotting = GRangesList(shortRange[[i]])
  for(j in 1:length(plotting)){
  plotting[[j]]$k = paste(plotting[[j]]$cluster,"binShort",j, sep = ":")
  }
  combinedPlotting[[i]] =  plotting
}

longRange
for(i in 1:length(longRange)){
  plotting = GRangesList(longRange[[i]])
  for(j in 1:length(plotting)){
    plotting[[j]]$k = paste(plotting[[j]]$cluster,"binLong",j, sep = ":")
  }
  combinedPlotting[[i]] =  c(combinedPlotting[[i]],plotting)
}





#save(combinedPlotting, file = "./combinedClustering.RData")





##############################
# Now plot the contact maps and 
#get coordinates of the clusters
##############################


matList = list()
clusterPositionsList = list()


for(j in 1:length(combinedPlotting)){
  
  
plotting =  unlist(combinedPlotting[[j]])

lengths = aggregate(mcols(plotting)$cluster, by = list(mcols(plotting)$k), FUN = length)
row.names(lengths) = lengths$Group.1
#for each cluster get the min start and max end
plottingSplit = split(plotting, paste(mcols(plotting)$k, mcols(plotting)$type))

minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
clusterPositionsList[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                              "ls" = minStarts[seq(1,length(minStarts),2)],
                              "le" = maxEnd[seq(1,length(maxEnd),2)],
                              "rs" = minStarts[seq(2,length(minStarts),2)],
                              "re" = maxEnd[seq(2,length(maxEnd),2)],
                              "size" = lengths[sub("\\s.*","",names(maxEnd)[seq(1,length(minStarts),2)]),])




matList[[j]] = matrix(0,nrow = 29903, ncol = 29903)


for(i in 1:nrow(clusterPositionsList[[j]])){
  matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
               clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] =    matList[[j]][clusterPositionsList[[j]][i,"ls"]:clusterPositionsList[[j]][i,"le"],
                                                                                                      clusterPositionsList[[j]][i,"rs"]:clusterPositionsList[[j]][i,"re"]] + clusterPositionsList[[j]][i, "size.x"]
  print(i)
}

# Un comment to save
#myCol = colorRampPalette(c("black",(brewer.pal(9,"YlOrRd"))))(30)
#ss <- ResizeMat(matList[[j]], c(10000,10000)) 
#pdf(paste("../clustering/clusterTestFullsars_S",j,".pdf"), height = 5, width = 5)
#heatmap3((log2((ss+1))), col=myCol, scale="none" ,Rowv = NA, Colv = NA)
#dev.off()
#write.table(clusterPositionsList[[j]],paste( "../clustering/SARS-CoV-2_cluster_positions_S",j,".txt"), sep = "\t",row.names = F, quote = F)

}





##############################
# Now cluster the clusters
##############################
clusterPositionsList
combinedPlotting

superclustersPoisitonList = list()
superclustersPlotting = list()


for(z in 1:4){
  
  clusterPositions = clusterPositionsList[[z]]
  # changes ranges that overlap
  clusterPositions2 = clusterPositions
  for(i in 1:nrow(clusterPositions )){
    if(clusterPositions$le[i] > clusterPositions$rs[i]){
      clusterPositions2[i,"rs"] =     clusterPositions[i,"le"] +1
    }
  }
  
  #make Granges
  left = GRanges(seqnames=IDSars,  
                 IRanges( 
                   start=clusterPositions2$ls, 
                   end=clusterPositions2$le 
                 )) 
  names(left) <- clusterPositions2$id 
  
  
  right= GRanges(seqnames=IDSars,  
                 IRanges( 
                   start=clusterPositions2$rs, 
                   end=clusterPositions2$re 
                 )) 
  names(right) <- clusterPositions2$id 
  
  
  
  
  #make a GRanges from the right 
  distances = GRanges(seqnames=IDSars,  
                      IRanges( 
                        start=clusterPositions2$le, 
                        end=clusterPositions2$rs 
                      )) 
  names(distances) <- clusterPositions2$id  
  
  
  # Now make super clusters
  adjacancyMat = getAdjacancyMat(distances,"nucleotide", 4)
  net = graph_from_adjacency_matrix(adjacancyMat, mode = "undirected", weighted = T)
  clustering = cluster_walktrap(net,steps = 1)
  highest_clusters = names(table(membership(clustering)))
  superclustersPlotting[[z]]  = printClustersFast("../clustering/combined/",clustering, highest_clusters, left, right)
  plottingListFull = superclustersPlotting[[z]]
  # Now for any cluster that didn't make it, add them to the file 
  length(unique(names(plottingListFull)))
  length(which(unique(names(plottingListFull)) %in% clusterPositions$id))
  missing = as.character(clusterPositions$id[which( !(as.character(clusterPositions$id) %in% unique(names(plottingListFull)) ) )])
  clusterPositionsmissing = clusterPositions[clusterPositions$id %in% missing,]
  superclusters = plottingListFull
  plottingListFull = superclusters 
  cluster = mcols(plottingListFull)$cluster
  names(cluster)= names(plottingListFull) 
  clusterPositions$superCluster = cluster[as.character((clusterPositions$id))]
  clusterPositions = clusterPositions[!is.na(clusterPositions$superCluster),]
  lengths = aggregate(clusterPositions$size.x, by = list(clusterPositions$superCluster), FUN = sum)
  row.names(lengths) = lengths$Group.1
  #for each cluster get the min start and max end
  plottingSplit = split(plottingListFull, paste(mcols(plottingListFull)$cluster, mcols(plottingListFull)$type))
  
  minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
  maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
  clusterPositionsCombined = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                        "ls" = minStarts[seq(1,length(minStarts),2)],
                                        "le" = maxEnd[seq(1,length(maxEnd),2)],
                                        "rs" = minStarts[seq(2,length(minStarts),2)],
                                        "re" = maxEnd[seq(2,length(maxEnd),2)],
                                        "size" = lengths[as.numeric(sub("\\s.*","",names(maxEnd)[seq(1,length(minStarts),2)])),])
  
  
  colnames(clusterPositionsmissing)
  colnames(clusterPositionsCombined)
  superclustersPoisitonList[[z]] = rbind.data.frame(clusterPositionsmissing,clusterPositionsCombined, stringsAsFactors = F)
  clusterPositions = superclustersPoisitonList[[z]]
  mat = matrix(0,nrow = 29903, ncol = 29903)
  for(i in 1:nrow(clusterPositions)){
    mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
        clusterPositions[i,"rs"]:clusterPositions[i,"re"]] =    mat[clusterPositions[i,"ls"]:clusterPositions[i,"le"],
                                                                    clusterPositions[i,"rs"]:clusterPositions[i,"re"]] + clusterPositions[i, "size.x"]
    print(i)
  }
  
  
  #uncomment to get heatmaps of clusters and tables
  #ss <- ResizeMat(mat, c(10000,10000)) 
  #myCol = colorRampPalette(c("black",(brewer.pal(9,"YlOrRd"))))(50)
  #pdf(paste("../clustering/superClusters",z,".pdf"), height = 6, width = 6)
  #heatmap3((log2((ss+1))), col=myCol, scale="none" ,Rowv = NA, Colv = NA)
  #dev.off()
  #write.table(superclustersPoisitonList[[z]], paste("../clustering/superClusterspositions",z,".txt"), sep = "\t",row.names = F, quote = F)
  
}

#un comment to save data
#save(clusterPositionsList,
#combinedPlotting,
#superclustersPoisitonList ,
#superclustersPlotting ,file = "clustersAndSuperClusters.rData")




# Trim the clusters 
combinedPlotting
clusterPositionsList
superclustersPlotting
superclustersPoisitonList

#first make a plotting GRanges list of each range in 
# each cluster 
allChimerasForSuperClustersPlotting = list()

for (i in 1:4){
  combinedPlottingSplit = combinedPlotting[[i]]
  combinedPlottingUnlist = unlist(combinedPlottingSplit)
  library(stringr)
  superClusterArray =  sub("\\s.*","",names(superclustersPlotting[[i]][duplicated(names(superclustersPlotting[[i]]))]))
  names(superClusterArray) = superclustersPlotting[[i]][duplicated(names(superclustersPlotting[[i]]))]$cluster
  x = superclustersPoisitonList[[i]][ grep("bin", row.names(superclustersPoisitonList[[i]])),]
  names = c(names(superClusterArray), row.names(x))
  superClusterArray = c(superClusterArray,row.names(x))
  names(superClusterArray) = names
  combinedPlottingUnlist$superCluster = "X"
  for( z in 1:length(superClusterArray)){
    supercluster = names(superClusterArray)[z]
    cluster = unique(sub("\\s.*","",superClusterArray[z]))
    combinedPlottingUnlist[combinedPlottingUnlist$k == cluster,]$superCluster = supercluster
  }
  allChimerasForSuperClustersPlotting[[i]] = combinedPlottingUnlist
}
allChimerasForSuperClustersPlotting
save(allChimerasForSuperClustersPlotting, superclustersPoisitonList,file =  "superClustersBeforeTrimming.Rdata")



##############################
# Now trim the clusters
##############################
i = 4
load("superClustersBeforeTrimming.Rdata")
allChimerasForSuperClustersPlottingTrimmed = list()
for(i in 1:4){
  allChimerasForSuperClustersPlottingTrimmed[[i]] = GRanges()
  for(cluster in unique(allChimerasForSuperClustersPlotting[[i]]$superCluster)){
    cluster2 = sub("\\s.*","" ,cluster)
    lefty = list()
    for(l in c("left","right")){
      clusterrange = allChimerasForSuperClustersPlotting[[i]][allChimerasForSuperClustersPlotting[[i]]$superCluster == cluster & allChimerasForSuperClustersPlotting[[i]]$type == l  ,]
      min = min(start(clusterrange[clusterrange$superCluster == cluster,]))
      max = max(end(clusterrange[clusterrange$superCluster == cluster,]))
      s = (start(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
      e = (end(clusterrange[clusterrange$superCluster == cluster &clusterrange$type == l  ,]))
      seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
      x = unlist(c(seq2(from = (s), to = e)))
      xt = table(x)
      x1 = x
      removal =  mean(x) + sd(x)*1
      included  = GRanges(seqnames=IDSars,  
                          IRanges( 
                            start=rep(min, length(clusterrange)), 
                            end=rep(removal, length(clusterrange))
                          )) 
      if( l == "left"){
        removal =  mean(x) -  sd(x)*1
        included = GRanges(seqnames=IDSars,  
                           IRanges( 
                             start=rep(removal, length(clusterrange)), 
                             end=rep(max, length(clusterrange))
                           )) 
      }

      t = pintersect( clusterrange,included)
      allChimerasForSuperClustersPlottingTrimmed[[i]] = c(  allChimerasForSuperClustersPlottingTrimmed[[i]],t)
      s = (start(t))
      e = (end(t))
      seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
      x = unlist(c(seq2(from = s, to = e)))
      if( l == "left"){
        lefty[[1]] = x
        lefty[[2]] = x1
      }
    }
    # un comment to print the views of the trimming
    #print(cluster )
    #tbl1 = data.frame(table(c(x1,lefty[[2]])))
    #tbl2 = data.frame(table(c(x,lefty[[1]])))
    #plot(ggplot(mapping =  aes(x = Var1, y = as.numeric(as.character(Freq))))+ 
    #       geom_bar(data = tbl1, stat = "identity")+
    #       geom_bar(data = tbl2, stat = "identity", colour = "firebrick") +
    #       theme_classic())
  }
}


allChimerasForSuperClustersPlotting
#save(allChimerasForSuperClustersPlottingTrimmed, file = "trimmedSuperCusters.Rdata")



# Now make a new plotting positions tables and plot
clusterPositionsListTrimmed = clusterPositionsList
matListTrimmed = matList
for(j in 1:4){


plotting  = allChimerasForSuperClustersPlottingTrimmed[[j]]

lengths = aggregate(mcols(plotting)$superCluster, by = list(mcols(plotting)$superCluster), FUN = length)
row.names(lengths) = lengths$Group.1


#for each cluster get the min start and max end
plottingSplit = split(plotting, paste(mcols(plotting)$superCluster, mcols(plotting)$type))

minStarts = unlist(lapply(plottingSplit, function(x) {return(min(start(x)))  }))
maxEnd = unlist(lapply(plottingSplit, function(x) {return(max(end(x)))  }))
x = sub("\\sleft","",names(maxEnd)[seq(1,length(minStarts),2)])
x = sub("\\sright","",x,2)
clusterPositionsListTrimmed[[j]] = data.frame("id" = names(maxEnd)[seq(1,length(minStarts),2)],
                                       "ls" = minStarts[seq(1,length(minStarts),2)],
                                       "le" = maxEnd[seq(1,length(maxEnd),2)],
                                       "rs" = minStarts[seq(2,length(minStarts),2)],
                                       "re" = maxEnd[seq(2,length(maxEnd),2)],
                                       "size" = lengths[x,])

matListTrimmed[[j]] = matrix(0,nrow = 29903, ncol = 29903)


for(i in 1:nrow(clusterPositionsListTrimmed[[j]])){

  matListTrimmed[[j]][clusterPositionsListTrimmed[[j]][i,"ls"]:clusterPositionsListTrimmed[[j]][i,"le"],
               clusterPositionsListTrimmed[[j]][i,"rs"]:clusterPositionsListTrimmed[[j]][i,"re"]] =    matListTrimmed[[j]][clusterPositionsListTrimmed[[j]][i,"ls"]:clusterPositionsListTrimmed[[j]][i,"le"],
                                                                                                      clusterPositionsListTrimmed[[j]][i,"rs"]:clusterPositionsListTrimmed[[j]][i,"re"]] + clusterPositionsListTrimmed[[j]][i, "size.x"]
  print(i)
  
}

# un comment to get heatmaps and tables
#myCol = colorRampPalette(c("black",(brewer.pal(9,"YlOrRd"))))(30)
#ss <- ResizeMat(matListTrimmed[[j]], c(10000,10000)) 
#pdf(paste("../clustering/clusterTestFullsars_TRIMMED_clusters_S",j,".pdf"), height = 5, width = 5)
#heatmap3((log2((ss+1))), col=myCol, scale="none" ,Rowv = NA, Colv = NA)
#dev.off()
#write.table(clusterPositionsListTrimmed[[j]],paste( "../clustering/SARS-CoV-2_cluster_positions_TRIMMED_clusters_S",j,".txt"), sep = "\t",row.names = F, quote = F)

}





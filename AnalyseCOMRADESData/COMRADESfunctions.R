#####################################################
#Functions for the analysis of COMRADES data 
#
#####################################################



 
rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){
  if(!require(fields)) stop("`fields` required.")
  
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  
  # interpolation
  ans[ncord] <- interp.surface(obj, loc)
  
  ans
}




#####################################################
# Load in the Files 
#####################################################

#load files from a directy cotaining the Hyb files

# takes a directory

#returns the list of hyb files 

readHybFiles = function(dir){
  list.files(dir, pattern = "hybrids.hyb$")
  #read in the files for analysis 
  hybList = list()
  
  
  
  for(i in 1:length(list.files(dir, pattern = "hybrids.hyb$"))){
    print(paste("Reading ",list.files(dir, pattern = "hybrids.hyb$")[i] ))
    sampleHyb = read.table(paste(dir,list.files(dir, pattern = "hybrids.hyb$")[i], sep =  ""), header = F, stringsAsFactors = F)
    hybList[[i]] = unique(sampleHyb)
    print(nrow(hybList[[i]]))
  }
  return(hybList)
}




#####################################################
# Subset HybList
#####################################################


subsetHybList = function(hybList, rna, number){
  RNA = rna
  
  smallest = number
  #test how many 18S sequences in each of the samples
  for(i in 1:length(hybList)){
    TE = RNA
    hybOutputO = hybList[[i]]
    hybOutput = hybOutputO[as.character(hybOutputO$V4) == rna & as.character(hybOutputO$V10) == rna,]
    nrowHyb = nrow(hybOutput )

    print(nrowHyb)
    
    sample = sample(1:nrowHyb, smallest, replace = F)
    hybList[[i]] = hybOutput[sample,]
    
  }
  

return(hybList)
}






#####################################################
# plot the matrices
#####################################################

plotMatrixList = function(hybMatList, directory, RNA, dir, a, b, c, d,h){ 
  require(RColorBrewer)
  require(heatmap3)
  for(hyb in 1:length(hybMatList[[1]])){
    
    
    for(TE in c(RNA)){
      
      print(list.files(dir, pattern = "hybrids.hyb$")[hyb])
      
      myCol = colorRampPalette(c("black",brewer.pal(9,"YlOrRd")))(50)

      #myCol <-  c("white",colorRampPalette(brewer.pal(11,"Spectral"))(1000))
      # myCol <- colorRampPalette(c("white", rep("black", 10000)))(10001)
      #if(hyb %%2 != 0){
      cols = log2(max(hybMatList[[TE]][[hyb]][a:b,c:d]+1)) / 12
      
      myCol = myCol[1:round(length(myCol) * cols, digits = 2)]
      #}
      print(myCol)
      print(cols)
      print(log2(max(hybMatList[[TE]][[hyb]][a:b,c:d]+1)) )
      print(round(length(myCol) * cols, digits = 2))
      pdf(paste(directory,list.files(dir, pattern = "hybrids.hyb$")[hyb] ,"_",TE,".pdf", sep = ""), height = h, width = h)
      heatmap3((log2(t(hybMatList[[TE]][[hyb]][a:b,c:d]+1))), col=myCol, scale="none" ,Rowv = NA, Colv = NA, useRaster = T)
      dev.off()
      
      
      
    }
  }
  

}






plotMatrixListfull = function(hybMatList, directory, RNA, dir,h){ 
  require(RColorBrewer)
  require(heatmap3)
  for(hyb in 1:length(hybMatList[[1]])){
    
    
    for(TE in c(RNA)){
      
      print(list.files(dir, pattern = "hybrids.hyb$")[hyb])
      
      myCol = colorRampPalette(c("black",brewer.pal(9,"YlOrRd")))(50)
      
      #myCol <-  c("white",colorRampPalette(brewer.pal(11,"Spectral"))(1000))
      # myCol <- colorRampPalette(c("white", rep("black", 10000)))(10001)
      if(hyb %%2 != 0){
        cols = log2(max(hybMatList[[TE]][[hyb]]+1)) / 25
        
        myCol = myCol[1:round(length(myCol) * cols, digits = 2)]
      }
      print(myCol)
      print(cols)
      print(log2(max(hybMatList[[TE]][[hyb]]+1)) )
      print(round(length(myCol) * cols, digits = 2))
      pdf(paste(directory,list.files(dir, pattern = "hybrids.hyb$")[hyb] ,"_",TE,".pdf", sep = ""), height = h, width = h)
      heatmap3((log2(t(hybMatList[[TE]][[hyb]]+1))), col=myCol, scale="none" ,Rowv = NA, Colv = NA, useRaster = T)
      dev.off()
      
      
      
    }
  }
  
  
}





#rescale and siz matrces

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){
  if(!require(fields)) stop("`fields` required.")
  
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  
  # interpolation
  ans[ncord] <- interp.surface(obj, loc)
  
  ans
}




#####################################################
# Make line plots
#####################################################


lineTraces = function(hybMatList, RNA, dir, rowOrCol){
  roundUp <- function(x) 10^trunc(log10(x))
  
  name = 1
  
  for(TE in c(RNA)){
    data2 = data.frame()
    for(hyb in 1:length(hybMatList[[TE]])){
      sampleName= hyb
      if(rowOrCol =="row"){      data = as.data.frame(rowSums(hybMatList[[TE]][[hyb]]))
      name = "row"
      }else if(rowOrCol =="col"){      data = as.data.frame(colSums(hybMatList[[TE]][[hyb]]))
      name = "col"}else if(rowOrCol =="both"){      
        
        d1 = as.data.frame(colSums(hybMatList[[TE]][[hyb]]))
        d2 = as.data.frame(rowSums(hybMatList[[TE]][[hyb]]))
        colnames(d1) = colnames(d2)
                           
        data = rbind.data.frame( d1,
                                 d2, stringsAsFactors = F )
      name = "both"
      } else{ print("please select either row or col")}
      data$pos = row.names(data)
      data$sample = hyb
      
      colnames(data) = c("chimeras","pos","sample")
      data2 = rbind(data,data2)
    }
    
    max = round(max(as.numeric(as.character(data2$pos))),-3)
    pdf(paste(dir,name ,"_",TE,".pdf", sep = ""), height = 2.5, width = 7.5)
    plot(ggplot()+
           geom_vline(xintercept = seq(0,max,1000), colour = "grey", alpha = 0.6)+
           geom_vline(xintercept = seq(0,max,500), colour = "grey", alpha = 0.4, linetype = "dashed")+
           geom_vline(xintercept = seq(0,max,100), colour = "grey", alpha = 0.3, linetype = "dotted")+
           geom_area(data = data2,aes(x = as.numeric(as.character(pos)), y = as.numeric(as.character(chimeras)), fill = sample))+
           facet_grid(sample~.)+
           theme_classic())+
      theme(legend.position = "none")
    
    dev.off()
    
    
    
  }

}












#####################################################
# Swap hyb files to make only those appearing on
# 3 prime first and also get only one RNA
#####################################################
#takes the hyb list produced by readHybFiles and an RNA

#outpute a new list with only that RNA and alll in one orientation

swapHybs = function(hybList, rna){
  
  for(hyb in 1:length(hybList)){
    
    #hybOutputO = hybList[[hyb]]
    hybList18S  =  hybList[[hyb]][hybList[[hyb]]$V4 == rna | hybList[[hyb]]$V10 == rna,]
    
    
    tmp1 = hybList18S[hybList18S$V7 < hybList18S$V13,]
    tmp2 = hybList18S[hybList18S$V7 > hybList18S$V13,]
    tmp2 = tmp2[,c("V1","V2","V3",
                   "V4","V11","V12",
                   "V13","V14","V15",
                   "V10","V5","V6",
                   "V7","V8","V9")]
    colnames(tmp2) = colnames(tmp1)
    comb = rbind.data.frame(tmp1, tmp2)
    
    
    
    
    hybList[[hyb]] = comb
  }
  
  return(hybList)
}


swapHybs2 = function(hybList, rna){
  
  for(hyb in 1:length(hybList)){
    
    hybList18S = hybList[[hyb]][as.character(hybList[[hyb]]$V4) == rna & as.character(hybList[[hyb]]$V10) == rna,]
    
    
    tmp1 = hybList18S[hybList18S$V7 < hybList18S$V13,]
    tmp2 = hybList18S[hybList18S$V7 > hybList18S$V13,]
    tmp2 = tmp2[,c("V1","V2","V3",
                   "V4","V11","V12",
                   "V13","V14","V15",
                   "V10","V5","V6",
                   "V7","V8","V9")]
    colnames(tmp2) = colnames(tmp1)
    comb = rbind.data.frame(tmp1, tmp2)
    
    
    
    
    hybList[[hyb]] = comb
  }
  
  return(hybList)
}



##############################
# Make the matrices
##############################

#takes a hyb list produced by either :
# readHybFiles or swapHybs along with the rna selected

#returns a list of the same dimensions containing the matrices of the chimeras 

getMatrices = function(hybList, rna){
  hybMatList = list()
  for(hyb in 1:length(hybList)){
    
    hybOutputO = hybList[[hyb]]
    print(hyb)
    for(TE in c(rna)){
      
      
      
      
      hybOutput =  hybOutputO[as.character(hybOutputO$V4) == rna & as.character(hybOutputO$V10) == rna,]
      hybOutput = unique(hybOutput)
      
      startsends = hybOutput[,c(7,8,13,14)]
      
      
      
      mat = matrix(0, ncol = max(startsends), nrow =max(startsends))
      for(i in 1:nrow(startsends)){
        data = startsends[i,]
        xData = seq(data$V7, data$V8)
        yData = seq(data$V14, data$V13)
        
        #mat[xData:yData,xData:yData] =         mat[xData:yData,xData:yData] +1
        
        for(j in 1:length(xData)){
          mat[xData[j],yData[j]] = mat[xData[j],yData[j]] +1
        }
      }
      
      hybMatList[[TE]][[hyb]] = mat
    }
  }
  return(hybMatList)
}










##############################
# Subset Hyb
##############################

# This function allows the dubsetting of the distance between 
# left and right sides of the chimeras

#takes a hyb list and returns a subsetted hyb list 

subsetHybList = function(hybList, min, max, length){
  longDistHyb = list()
  for (i in 1:length(hybList)){
    hybList[[i]]$dist = hybList[[i]]$V13 -  hybList[[i]]$V8
    longDistHyb[[i]] = hybList[[i]][hybList[[i]]$dist < max & hybList[[i]]$dist >= min,]
    leftLength = longDistHyb[[i]]$V6 - longDistHyb[[i]]$V5
    rightLength = longDistHyb[[i]]$V12 - longDistHyb[[i]]$V11
    longDistHyb[[i]] = longDistHyb[[i]][leftLength < length & rightLength < length,]
  }
  return(longDistHyb)
}







##############################
# Make the Granges
##############################

#takes a hyb format files and returns GrangesList of the left right and the gap



hybToGRanges = function(hybOutput2, seqName){
  
  gList = list()
  for(i in 1:length(hybOutput2)){
    hybOutput = hybOutput2[[i]]
    gList[[i]] = GRangesList()
    #make a Granges from the left
    left <- GRanges(seqnames=seqName, 
                    IRanges( 
                      start=hybOutput$V7, 
                      end=hybOutput$V8
                    )) 
    names(left) <- hybOutput$V1 
    
    
    
    
    #make a GRanges from the right 
    right <-GRanges(seqnames=seqName, 
                    IRanges( 
                      start=hybOutput$V13, 
                      end=hybOutput$V14 
                    )) 
    names(right) <- hybOutput$V1 
    
    
    
    distances = GRanges(seqnames=seqName,  
                        IRanges( 
                          start=end(left), 
                          end=start(right) 
                        )) 
    names(distances) <- hybOutput$V1 
    
    gList[[i]][["left"]] = left
    gList[[i]][["right"]] = right
    gList[[i]][["gap"]] = distances
  }
  return(gList)
}











##############################
# Make and adjancy matrix of theoverlaps 
##############################
#Make a weights adjacancy matrix

 
getAdjacancyMat = function(hybGranges, nucletideOrPerc, cutoff){
  distances = hybGranges
  max = max(width(distances))
  
  #get overlapping
  hits <- findOverlaps(distances, drop.self=T, drop.redundant=F)
  # get the relative overlap for the weights 
  x <- distances[queryHits(hits)]
  y <- distances[subjectHits(hits)]
  
  relative_overlap <-  width(pintersect(x, y)) / pmax(width(x), width(y))  

  hitsWithOverlap = hits
  # parameter for relative overlap 
  
  print(length(hitsWithOverlap))
  
  if(nucletideOrPerc == "none"){
    hitsWithOverlap = hits
  } else if(nucletideOrPerc == "nucleotide"){
    relative_overlap =     pmax(width(x), width(y)) - width(pintersect(x, y))
    relative_overlap = cutoff - relative_overlap
    
    #relative_overlap =  width(pintersect(x, y)) 
    hitsWithOverlap = hits[relative_overlap <= cutoff & relative_overlap >= 0 ]
    relative_overlap = relative_overlap[relative_overlap <= cutoff  & relative_overlap >= 0 ]
    #relative_overlap = relative_overlap / pmax(width(x), width(y))
    #relative_overlap = relative_overlap2
    #relative_overlap2 =  sqrt((relative_overlap -(max(relative_overlap)+1))^2)
  } else if(nucletideOrPerc == "perc"){
    #relative_overlap <-  width(pintersect(x, y)) / pmax(width(x), width(y))  
    relative_overlap = (1- (width(pintersect(x, y)) / max))
    #points(pmax(width(x), width(y)), relative_overlap , col =width(pintersect(x, y)) )
    #relative_overlap = relative_overlap / as.numeric(paste("1.", (pmax(width(x), width(y))), sep = "")) 
    hitsWithOverlap = hits[relative_overlap >= cutoff]
    if(length(hitsWithOverlap ) ==0 ){return(0)}else{}
    print(length(hitsWithOverlap))
    relative_overlap = relative_overlap[relative_overlap >= cutoff]

  }
  
  hitsMat = as.data.frame(hitsWithOverlap)
  hitsMat$weight = relative_overlap 
  #hitsMat2 = hitsMat[,c("subjectHits","queryHits","weight")]
  #colnames(hitsMat2) = colnames(hitsMat)
  #hitsMat = rbind.data.frame(hitsMat,hitsMat2)

  testLong = dcast(hitsMat, queryHits ~ subjectHits, value.var = "weight")
  row.names(testLong) = testLong$queryHits
  testLong = testLong[,-1]
  testLong = as.matrix(testLong)
  testLong[is.na(testLong)] =0
  
  return(testLong)

}



##############################
# Plot Chimeras
##############################

plotChimeras = function(chimeraList, subset){
# combine the granges
  combined = c(chimeraList[[1]]$left[1:subset],chimeraList[[1]]$right[1:subset])
  combined 
  # plt them 
  options(ucscChromosomeNames=FALSE) 
  
  
  gtrack <- GenomeAxisTrack(add53=T,labelPos="above",grid=T)
  annoT <- AnnotationTrack(combined, 
                           group = c(names(chimeraList[[1]]$left[1:subset]),names(chimeraList[[1]]$left[1:subset]))) 
  
  
  plotTracks(list(gtrack,annoT)) 

}







##############################
# Plot the clusters and create a granges of the clusters 
##############################
#needs as input dir, clustering, highest_clusters, left, right
# then output a number of plots in the directory

# plus the granges of the information for the hghest clusters. 


#helper function for below
addCluster = function(granges, indexes, prev, cluster, type){
  x = granges[as.numeric(indexes)]
  x$cluster = cluster
  x$type = type
  prev= c(prev, x)
}




printClusters = function(dir, clustering, highest_clusters, left, right ){
  
  plotting = GRanges()
  for(i in highest_clusters ){#:(max(as.numeric((names(table(membership(cluster3))))))-300)){
    c1c1 = names(membership(clustering)[membership(clustering) == i])
    #c1c1 = which(clustering == i)
    
    #pdf("testMatric.pdf", height = 20, width = 20)
    #heatmap3(subset, scale="none",showColDendro = T,showRowDendro = T)
    #dev.off()
    plotting2 = GRanges()
    plotting2 = addCluster(left, c1c1, plotting2,i  ,"left")
    plotting2 = addCluster(right, c1c1, plotting2,i ,"right" )
    annoT <- AnnotationTrack(plotting2,group = names(plotting2) , fill = plotting2$cluster, stacking = "squish",chromosome = "NC045512_hsapiens_SARSCOV-1_SARSCOV", col=NULL) 
    gtrack <- GenomeAxisTrack(add53=T,labelPos="above",grid=T)
    #from = min(start(plotting2))+1000
    #to = max(start(plotting2)) - 1000
    pdf(paste(dir,"cluster_",i,"_",min(start(plotting2)),".pdf", sep = ""), height = (length(plotting2)/50)+3)
    plotTracks( list(gtrack,annoT), groupingAnnotation="group", chromosome = "NC045512_hsapiens_SARSCOV-1_SARSCOV") 
    dev.off()
    
    #ad the info to the GRanges
    # plotting = addCluster(distances, c1c1, plotting,i,"middle"  )
    plotting = addCluster(left, c1c1, plotting,i  ,"left")
    plotting = addCluster(right, c1c1, plotting,i ,"right" )
    
  }
#   gtrack <- GenomeAxisTrack(genome,chromosme = "18S")
  annoT <- AnnotationTrack(plotting,group = names(plotting),stacking = "squish", fill = plotting$cluster,chromosome = "NC045512_hsapiens_SARSCOV-1_SARSCOV", col=NULL,id=plotting$cluster) 
  
   pdf(paste(dir,"cluster_all.pdf", sep = ""), height = (length(plotting)/58)+2)
   gtrack <- GenomeAxisTrack(add53=T,labelPos="above",grid=T)
  plotTracks( list(gtrack,annoT), groupingAnnotation="group", chromosome = "NC045512_hsapiens_SARSCOV-1_SARSCOV", showFeatureId=T,cex = 0.1) 
  dev.off()
  return(plotting)
}










printClustersFast = function(dir, clustering, highest_clusters, left, right ){
  
  plotting = GRanges()
  for(i in highest_clusters ){#:(max(as.numeric((names(table(membership(cluster3))))))-300)){
    c1c1 = names(membership(clustering)[membership(clustering) == i])
    plotting2 = GRanges()
    plotting2 = addCluster(left, c1c1, plotting2,i  ,"left")
    plotting2 = addCluster(right, c1c1, plotting2,i ,"right" )
    plotting = addCluster(left, c1c1, plotting,i  ,"left")
    plotting = addCluster(right, c1c1, plotting,i ,"right" )
    
  }

  return(plotting)
}



#get gaps from clusters

getGapsFromClusters = function(clustering, highest_clusters, gaps ){
  
  plotting = GRanges()
  for(i in highest_clusters ){#:(max(as.numeric((names(table(membership(cluster3))))))-300)){
    c1c1 = names(membership(clustering)[membership(clustering) == i])

    plotting = addCluster(gaps, c1c1, plotting,i ,"gap" )
    
  }

  return(plotting)
}













##############################
# Plots the istribution of the clusters
##############################


plotClusterDist = function(clusters, known, plotting){ 
  boths = data.frame()
  for(i in 1:length(clusters)){
    cluster = clusters[i]
    x = melt(as.vector(coverage(plotting[plotting$cluster == cluster & plotting$type == "right"])[[1]]))
    x$x = row.names(x)
    x$side = "right"
    y = melt(as.vector(coverage(plotting[plotting$cluster == cluster & plotting$type == "left"])[[1]]))
    y$x = row.names(y)
    y$side = "left"
    both = rbind.data.frame(x,y)
    both$cluster = cluster
    boths = rbind.data.frame(boths,both)
  }
  boths = boths[boths$value > 0, ]
  
  if(known == T){
    boths2 = data.frame()
    known = read.table("/Users/jp/projects/omerRiboHighTemp/newPipeline/riboVision_longDistanceKnown.txt", header = T)
    for(i in unique(known$cluster)){
      
      left <- GRanges(seqnames="18S", 
                      IRanges( 
                        start=known$lS, 
                        end=known$lE
                      ), cluster = known$cluster, side = "left") 
      
      #make a GRanges from the right 
      right <-GRanges(seqnames="18S", 
                      IRanges( 
                        start=known$Rs, 
                        end=known$Re 
                      ), cluster = known$cluster, side = "right") 
      x = melt(as.vector(coverage(left[left$cluster == i])[[1]]))
      x$x = row.names(x)
      x$side = "left"
      y = melt(as.vector(coverage(right[right$cluster == i])[[1]]))
      y$x = row.names(y)
      y$side = "right"
      both = rbind.data.frame(x,y)
      both$cluster =i
      both[both$value == 2,"value"] = 1
      boths2 = rbind.data.frame(boths2,both)
    }
    
    
    boths = rbind.data.frame(boths,boths2)
  }
  ggplot(data = boths, aes(x = as.numeric(x),y = value, colour = side)) +
    geom_vline(xintercept = seq(0,1850,10), alpha = 0.2, colour = "grey", linetype = "dotted")+ 
    geom_vline(xintercept = seq(0,1850,100), alpha = 0.4, colour = "grey")+ 
    geom_step()+
    theme_classic()+
  #  xlim(1650,1900)+
    facet_grid(cluster~., scales = "free_y")+
    xlim(0,1870)
}















##########################################################################################
#
# function to add known to the plot
#
##########################################################################################
know18Sints = read.table("ribovision18S.txt", header = F)

addKnown = function(known, table){
  known$index1 = paste(known[,1],known[,2], sep=":")
  known$index2 = paste(known[,2],known[,1], sep=":")
  
  table$index1 = paste(table[,1], table[,2], sep=":")
  
  
  table$isKnown = 0
  
  for(i in 1:nrow(table)){
    if(table$index1[i] %in% c(known$index1, known$index2)){
      table$isKnown[i] = 1
    }
  }
  return(table)
  
}



##########################################################################################
#
# function with sequences taken from the distribution at random lengths
#
##########################################################################################

getClusterSeqsFromDist = function(cluster){
  
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = plotting[plotting$cluster == cluster & plotting$type == i,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    covVect = covVect[covVect > (mean(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    
    seq = seq(min,max)
    mean = sum(seq * covVect) / sum(covVect)
    sigma = sqrt(sum((seq  - mean)**2 * covVect) / (sum(covVect)-1))
    coords[[i]] = c(min,max,mean,sigma)
    
  }
  
  
  
  
  
  startPos = list()
  endPos = list()
  seqs = list()
  for(i in c("left","right")){
    #grangesClusterfinal = plotting[plotting$cluster == cluster & plotting$type == i & start(plotting) > coords[[i]][1] & end(plotting) <  coords[[i]][2] ,]
    mids =   round(rnorm(20,coords[[i]][3],coords[[i]][4]))
    lengths = round(runif(100,5,20))
    starts = mids-lengths
    ends = mids+lengths
    for(j in 1:length(starts)){
      seqs[[i]][j] = paste(fasta$RN18S1[starts[j]:ends[j]], collapse = "")
    }
    startPos[[i]] = starts
    endPos[[i]] = ends
    
  }
  
  minSeqs = min(length(startPos[[1]]), length(startPos[[2]]))
  indexes = round(runif(40,1, minSeqs))
  
  return(list(startPos,endPos, indexes, seqs))
}





##########################################################################################
#
# function with random subset
#
##########################################################################################

getClusterClusterBinding = function(cluster, plotting){
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = plotting[plotting$cluster == cluster & plotting$type == i,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    #covVect = covVect[covVect > (mean(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    coords[[i]] = c(min,max)
    
  }
  
  #seq = paste(fasta$RN18S1[min:max], collapse = "")
  #command = paste("echo \">right\n",seq,"\" | RNAfold -p --MEA ", sep = "")
  #system(command)
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  #get the sequences from the most representative of the clusters
  startPos = list()
  endPos = list()
  seqs = list()
  for(i in c("left","right")){
    grangesClusterfinal = plotting[plotting$cluster == cluster & plotting$type == i & start(plotting) >= coords[[i]][1] & end(plotting) <=  coords[[i]][2] ,]
    starts = start(grangesClusterfinal)
    ends = end(grangesClusterfinal)
    for(j in 1:length(starts)){
      seqs[[i]][j] = paste(fasta$RN18S1[starts[j]:ends[j]], collapse = "")
    }
    startPos[[i]] = starts
    endPos[[i]] = ends
    
  }
  
  minSeqs = min(length(startPos[[1]]), length(startPos[[2]]))
  indexes = round(runif(10,1, minSeqs))
  
  return(list(startPos,endPos, indexes, seqs))
}


##########################################################################################
#
#Function with tilling
#
##########################################################################################


#now fold them

getClusterClusterBindingTilling = function(cluster){
  
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = plotting[plotting$cluster == cluster & plotting$type == i,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    covVect = covVect[covVect > (mean(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    coords[[i]] = c(min,max)
    
  }
  
  
  startPos = list()
  endPos = list()
  seqs = list()
  for(i in c("left","right")){
    #grangesClusterfinal = plotting[plotting$cluster == cluster & plotting$type == i & start(plotting) > coords[[i]][1] & end(plotting) <  coords[[i]][2] ,]
    starts = seq(coords[[i]][1],coords[[i]][2]-23,2)
    ends = starts+25
    for(j in 1:length(starts)){
      seqs[[i]][j] = paste(fasta$RN18S1[starts[j]:ends[j]], collapse = "")
    }
    startPos[[i]] = starts
    endPos[[i]] = ends
    
  }
  
  minSeqs = min(length(startPos[[1]]), length(startPos[[2]]))
  indexes = round(runif(50,1, minSeqs))
  
  return(list(startPos,endPos, indexes, seqs))
}



##########################################################################################
#
#   Function for short range with  a smal section from each side 
#
##########################################################################################


#now fold them
getClusterClusterShortRange= function(cluster, plotting, length){
  
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = plotting[plotting$cluster == cluster & plotting$type == i,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    covVect = covVect[covVect > (mean(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    pos = 0
    if(i == "left"){
      pos = max
    }else if(i == "right"){
      pos = min
    }
    
    coords[[i]] = pos
    
  }
  
  
  startPos = list()
  endPos = list()
  seqs = list()
  for(i in c("left","right")){
    #grangesClusterfinal = plotting[plotting$cluster == cluster & plotting$type == i & start(plotting) > coords[[i]][1] & end(plotting) <  coords[[i]][2] ,]
  
    starts = 0
    ends = 0
    if(i == "left"){
      ends = coords[[i]]
      starts = coords[[i]] - length
    }else if(i == "right"){
      starts = coords[[i]]
      ends = coords[[i]] + length
    }
    
    
    for(j in 1:length(starts)){
      seqs[[i]][j] = paste(fasta$RN18S1[starts:ends], collapse = "")
    }
    startPos[[i]] = starts
    endPos[[i]] = ends
    
  }
  
  minSeqs = min(length(startPos[[1]]), length(startPos[[2]]))
  indexes = 1
  
  return(list(startPos,endPos, indexes, seqs))
}







##########################################################################################
#
#   Function full sequence from each side 
#
##########################################################################################


#now fold them
getClusterClusterShortRangeWhole = function(cluster, plotting, length){
  
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = plotting[plotting$cluster == cluster,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    covVect = covVect[covVect > (mean(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    pos = 0
    if(i == "left"){
      pos = min
    }else if(i == "right"){
      pos = max
    }
    
    coords[[i]] = pos
    
  }
  
  
  seqs = paste(fasta$RN18S1[coords[[1]]:coords[[2]]], collapse = "")

  indexes = 1
  
  return(list(coords[[1]],coords[[2]], indexes, seqs))
}

























#now fold them
getClusterClusterShortRangeMidGap = function(cluster, clusterGap){
  
  library(seqinr)
  fasta = read.fasta("18s.fasta")
  
  coords = list()
  for(i in c("left","right")){
    grangesClusterRight = clusterGap[clusterGap$cluster == cluster,]
    
    
    covVect = as.vector(coverage(grangesClusterRight)[[1]])[as.vector(coverage(grangesClusterRight)[[1]]) !=0]
    covVect2 = as.vector(coverage(grangesClusterRight)[[1]])
    covVect = covVect[covVect > (mean(covVect)+sd(covVect))]
    min = min(which(covVect2 == covVect[1]))
    max = max(which(covVect2 == covVect[length(covVect)]))
    pos = 0
    if(i == "left"){
      pos = max
    }else if(i == "right"){
      pos = min
    }
    
    coords[[i]] = pos
    
  }
  
  
  startPos = list()
  endPos = list()
  seqs = list()
  for(i in c("left","right")){
    #grangesClusterfinal = plotting[plotting$cluster == cluster & plotting$type == i & start(plotting) > coords[[i]][1] & end(plotting) <  coords[[i]][2] ,]
    
    starts = 0
    ends = 0
    if(i == "left"){
      ends = coords[[i]]
      starts = coords[[i]] - 20
    }else if(i == "right"){
      starts = coords[[i]]
      ends = coords[[i]] + 20
    }
    
    
    for(j in 1:length(starts)){
      seqs[[i]][j] = paste(fasta$RN18S1[starts:ends], collapse = "")
    }
    startPos[[i]] = starts
    endPos[[i]] = ends
    
  }
  
  minSeqs = min(length(startPos[[1]]), length(startPos[[2]]))
  indexes = 1
  
  return(list(startPos,endPos, indexes, seqs))
}









##########################################################################################
#
# function to plot the bpbp interactions
# needs the start positions of startpos list left and right same for endpos alogn with 
# list of sequences left and right and the indexes of thos chosen
##########################################################################################



findBasePairsRNAfold = function(startPos, endPos, seqs,fasta){
  
  library(R4RNA)
      table = data.frame(stringsAsFactors = FALSE)
    #  for(i in 1:10){
      #command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold --ImFeelingLucky", sep = "")
      command = paste("echo \">",startPos,"-",endPos,"\n",seqs,"\" | RNAfold ", sep = "")
      x = system(command,intern = T)
      
      #extract the vienna and the start locations
      vienna = sub(" .*","",x[3])

      helix = viennaToHelix(vienna )
      
      helix = helix[,-c(3,4)]
      helix$rl = helix$i + startPos -1 
      helix$rr = helix$j + startPos -1 
      helix$pl = fasta[[1]][helix$i]
      helix$pr = fasta[[1]][helix$j]
  colnames(table) = c("l","r","rl","rr","pl","pr")
  table$check = 1
  aggTable = aggregate(table$check, by=list(table$rl,table$rr,table$pl,table$pr), FUN = sum )
  return(aggTable)
}























##########################################################################################
#
# function to plot the bpbp interactions
# needs the start positions of startpos list left and right same for endpos alogn with 
# list of sequences left and right and the indexes of thos chosen
##########################################################################################

startPos = list()
startPos[[1]] = startPos1
startPos[[2]] = startPos2
endPos = list()
endPos[[1]] = endPos1
endPos[[2]] = endPos2
seqs = list()
seqs[[1]] = seq1
seqs[[2]] = seq2

findBasePairs = function(startPos, endPos, indexes, seqs){
  
  library(R4RNA)
  table = data.frame(stringsAsFactors = FALSE)
  for(leftseq in indexes){
    for(rightseq in indexes){
      #command = paste("echo \">",startPos[[2]][i],"-",endPos[[2]][[i]],"-",startPos[[1]][[i]],"-",endPos[[1]][[i]],"\n",seqs[[2]][[i]],"&",seqs[[1]][[i]],"\" | RNAcofold -e 100 ", sep = "")
      #command = paste("echo \">",startPos[[1]][leftseq],"-",endPos[[1]][leftseq],"\n",seqs[[1]][leftseq],"\n>",startPos[[2]][rightseq],"-",endPos[[2]][rightseq],"\n",seqs[[2]][rightseq],"\" | RNAduplex -e 100", sep = "")
      
      #create vienna command and run capturing output
      command = paste("echo \">",startPos[[1]],"-",endPos[[1]],"\n",seqs[[1]],"\n>",startPos[[2]],"-",endPos[[2]],"\n",seqs[[2]],"\" | RNAduplex", sep = "")
      x = system(command,intern = T)
      
      #extract the vienna and the start locations
      vienna = sub(" .*","",x[3])
      startl = as.numeric(gsub(".* (\\d{1,2}),\\d{1,2}.*:.*(\\d{1,2}),\\d{1,2} .*", "\\1", x[3], perl = T))
      startr = as.numeric(gsub(".* (\\d{1,2}),\\d{1,2}.*:.*\\d{1,2},(\\d{1,2}) .*", "\\2", x[3], perl = T))
      #print(nchar(vienna))
      
      viennar =   startPos[[2]][rightseq] + (startr-1)
      viennal =   startPos[[1]][leftseq] + (startl-1)
      
      additionl = 0
      additionr = 0
      for(i in 1:nchar(vienna)){
        j = i+ additionl
        k = i + additionr 
        l = substr(vienna,j,j)
        r = substr(vienna,(nchar(vienna)+1)-k,(nchar(vienna)+1)-k)
        #print(paste(j,k,l,r))
        
        
        if(l == "(" & r == ")"){
          posl = (j-1) + viennal
          posr =  viennar - (k-1)
          p1 = fasta$RN18S1[posl]
          p2 = fasta$RN18S1[posr]
          #print(p1)
          #table = rbind.data.frame(table,c(j,k,posl,posr,fasta$RN18S1[posl],fasta$RN18S1[posr]))
          table = rbind.data.frame(table,c(as.numeric(j),as.numeric(k),as.numeric(posl),as.numeric(posr),as.character(p1),as.character(p2)),stringsAsFactors = FALSE)
        }
        
        if(l == "." & r == ")"){
          additionr = additionr -1
        }
        if(l == "(" & r == "."){
          additionl = additionl -1
        }
        if(l == "&" || r == "&"){
          break
        }
      }
    }
  }
  
  colnames(table) = c("l","r","rl","rr","pl","pr")
  table$check = 1
  aggTable = aggregate(table$check, by=list(table$rl,table$rr,table$pl,table$pr), FUN = sum )
  return(aggTable)
}












##########################################################################################
# Function that finds the ammount of evidence for a given interaction
##########################################################################################


getViennaForRegion = function(table, tableAll, output){
  
  minPos = min(c(table$i,table$j))
  minPos = minPos 
  maxPos = max(c(table$i,table$j))
  table$length = 1
  table$value = 1
  
  # table$i = table$i - minPos+1
  #  table$j = table$j - minPos+1
  #table
  #table = as.helix(table, length = (maxPos - minPos -1))
  
  table2 = table
  table2$i = table2$i - minPos+1
  table2$j = table2$j - minPos+1
  table2 = as.helix(table2, length = (maxPos - minPos -1))
  fseVienna = helixToVienna(table2)
  
  print(fseVienna)
  
  
  tableAllSample = tableAll[tableAll$sample %in% c("sample1","sample2"),]
  tableAllsarsS = tableAllSample
  tableAllsarsAgg = aggregate(tableAllsarsS$evidence, by = list(as.numeric(tableAllsarsS$p1), as.numeric(tableAllsarsS$p2)), FUN = sum)
  tableAllsarsAgg[order(tableAllsarsAgg$Group.1),]
  
  
  vec = tableAllsarsAgg$x
  names(vec) = paste(tableAllsarsAgg$Group.1, tableAllsarsAgg$Group.2, sep = ":")
  table$index = paste(table$i, table$j, sep = ":")
  table$evidence = vec[table$index]
  table$index2 = paste(table$j, table$i, sep = ":")
  table$evidence2 = vec[table$index]
  table[is.na(table$evidence2),'evidence2'] = 0
  
  x = rbind.data.frame(table[,c(1,8)],
                       setNames(table[,c(2,8)], names(table[,c(1,8)])),stringsAsFactors = F)
  x$evidence2 = log(x$evidence2+1)
  x = x[order(x$i),]
  x
  
  
  vec = x$evidence2
  names(vec) = x$i
  
  df = data.frame(minPos:maxPos,rep(0,length(minPos:maxPos)))
  colnames(df) = c("V1","V2")
  for(i in 1:nrow(df)){
    df$V2[i] = vec[as.character(df$V1[i])]
  }
  df[is.na(df$V2),2] = 0
  df$V1 = df$V1 - (min(df$V1)-1)
  df[is.na(df$V2),2] = 0
  write.table(df, output,quote = F, row.names = F, col.names = F)
  
  
}







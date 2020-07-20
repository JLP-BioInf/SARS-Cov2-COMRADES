#######################################################################################
#######################################################################################
#                                                                                     #
#    These Functions analyse COMRADES data - Author - Jonathan Price                  #
#                                                                                     #
#######################################################################################
#######################################################################################



#' Helper for Resize a Matrix 
#'
#' helper for Resize a Matrix to make it smaller for plotting
#' @param x 
#' @param newrange
#' @examples 
#' loc[,1] = rescale(ncord[,1], c(1,odim[1]))
#' @export
 
rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}


#' Resize a Matrix 
#'
#' Resize a Matrix to make it smaller for plotting
#' @param mat The matrix to be re-sized
#' @param ndim The dimensions to re-size to c(10,10)
#' @return The resized matrix
#' @examples 
#' ResizedMat <- ResizeMat(mat, c(10000),c(10000));
#' @export


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




#' Read in hyb files
#'
#' Read in hyb files and returns a list of hyb data
#' @param dir the directory that contains the *hybrids.hyb files
#' @return A list with the hyb files loaded 
#' @examples 
#' @authors Jonathan Price
#' hybList <- readHybFiles("/path/to/hyb/files/");
#' @export


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




#' Subset a list of hyb files
#'
#' Function used to subset a list of hyb data created by readHybFiles
#' This function produces the same size list as before but with an 
#' extra layer showing with the rna of interest. 
#' @param hybList the original hybList created with readHybFiles
#' @param rna the rna of interest that you want to subset 
#' @param number The number of randomly subsetted chimeric reads you need
#' @return A list of subsetted hyb files
#' @authors Jonathan Price
#' @examples 
#' hybListSubset <- subsetHybList(hybList, "myRNA" ,  10000);
#' @export



subsetHybList = function(hybList, rna, number){
  RNA = rna
  
  smallest = number
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









#' Make a matric of contacts interactions
#'
#' Function used to create a list of matrices for plotting with
#' plotMatrixList or plotMatrixListFull, the output list will be same as the 
#' input exect for an extra list layer for the specific RNA
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna the rna of interest that you want to subset 
#' @return A list of matrices
#' @authors Jonathan Price
#' @examples 
#' hybMatList <- subsetHybList(hybList, "myRNA");
#' @export


getMatrices = function(hybList, rna){
  hybMatList = list()
  for(hyb in 1:length(hybList)){
    
    hybOutputO = hybList[[hyb]]
    for(TE in c(rna)){
      
      
      hybOutput =  hybOutputO[as.character(hybOutputO$V4) == rna & as.character(hybOutputO$V10) == rna,]
      hybOutput = unique(hybOutput)
      
      startsends = hybOutput[,c(7,8,13,14)]
      
      
      
      mat = matrix(0, ncol = max(startsends), nrow =max(startsends))
      for(i in 1:nrow(startsends)){
        data = startsends[i,]
        xData = seq(data$V7, data$V8)
        yData = seq(data$V14, data$V13)
      
        
        for(j in 1:length(xData)){
          mat[xData[j],yData[j]] = mat[xData[j],yData[j]] +1
        }
      }
      
      hybMatList[[TE]][[hyb]] = mat
    }
  }
  return(hybMatList)
}







#' Plot a subsetted matrix List 
#'
#' Function used to plot a subsetted matrix of hyb datas created 
#' with getMatrices. Mainly used if RNA genome of interest is too large to 
#' create a usable heatmap
#' @param hybMatList list of contact matrices created with getMatrices
#' @param directory The path to the derectory where heatmaps should be printed
#' @param RNA The rna of interest
#' @param dir the directory that contains the *hybrids.hyb files
#' @param a x1 - coordinate
#' @param b x2 - coordinate
#' @param c y1 - coordinate
#' @param d y2 - coordinate
#' @param h the size of the plot in inches
#' @authors Jonathan Price
#' @examples 
#' plotMatrixList(hybMatList,"/path/to/output/, "myRNA","/path/to/hyb/files",10,10,10,10,5);
#' @export



plotMatrixList = function(hybMatList, directory, RNA, dir, a, b, c, d,h){ 
  require(RColorBrewer)
  require(heatmap3)
  for(hyb in 1:length(hybMatList[[1]])){
    for(TE in c(RNA)){
      
      print(list.files(dir, pattern = "hybrids.hyb$")[hyb])
      
      myCol = colorRampPalette(c("black",brewer.pal(9,"YlOrRd")))(50)

      cols = log2(max(hybMatList[[TE]][[hyb]][a:b,c:d]+1)) / 12
      
      myCol = myCol[1:round(length(myCol) * cols, digits = 2)]
      
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




#' Plot a subsetted matrix List 
#'
#' Function used to plot a subsetted matrix of hyb datas created 
#' with getMatrices. Mainly used if RNA genome of interest is too large to 
#' create a usable heatmap
#' @param hybMatList list of contact matrices created with getMatrices
#' @param directory The path to the derectory where heatmaps should be printed
#' @param RNA The rna of interest
#' @param dir the directory that contains the *hybrids.hyb files
#' @param h the size of the plot in inches
#' @authors Jonathan Price
#' @examples 
#' plotMatrixList(hybMatList,"/path/to/output/, "myRNA","/path/to/hyb/files/",5);
#' @export




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








#' Make a line plot of the colSums of a hybMatList
#'
#' Used to make a line plot of the colSums of a hybMatList, hybMatList can
#' be made with getMatrices
#' @param hybMatList list of contact matrices created with getMatrices
#' @param RNA The rna of interest
#' @param dir the directory where line plots should be printed
#' @authors Jonathan Price
#' @examples 
#' lineTraces(hybMatList,"myRNA", "/path/to/output/","row");
#' @export




lineTraces = function(hybMatList, RNA, dir, rowOrCol){
  roundUp <- function(x) 10^trunc(log10(x))
  
  name = 1
  
  for(TE in c(RNA)){
    data2 = data.frame()
    for(hyb in 1:length(hybMatList[[TE]])){
      sampleName= hyb
      if(rowOrCol =="row"){      
        data = as.data.frame(rowSums(hybMatList[[TE]][[hyb]]))
        name = "row"
      }else if( rowOrCol =="col" ){     
        data = as.data.frame(colSums(hybMatList[[TE]][[hyb]]))
        name = "col"
      } else { 
        print("please select either row or col")
      }
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











#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Also subsets to make sure that at least 
#' one side of the duplex is on your RNA of interest 
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of "swapped" hyb datas
#' @authors Jonathan Price
#' @examples 
#'  hybListSwapped =  swapHybs(hybList,"myRNA");
#' @export



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





#' Swap the table to ensure that 3 prime most duplex side is on the left of the table
#'
#' Swap the table to ensure that 3 prime most duplex side is ont he left of the table
#' used to make one sides heatmaps and other reasons where having the left of the table
#' coming after the right side is a problem. Different from swapHybs as it 
#' ensure that BOTH duplex sides originate from the RNA of interest.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of "swapped" hyb datas
#' @authors Jonathan Price
#' @examples 
#'  hybListSwapped =  swapHybs2(hybList,"myRNA");
#' @export


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








#' Makes Granges of HybLists
#'
#' This function is useful to turn a list of hyb datas into lists of GRanges
#' It creates a list for each sample one for the left side one for the right 
#' side and one for the gap in the middle.
#' @param hybList the original hybList created with readHybFiles or subsetHybList
#' @param rna The rna of interest
#' @return A list of GRanges data in hyb format
#' @authors Jonathan Price
#' @examples 
#'  hybGranges =  hybToGRanges(hybList,"myRNA");
#' @export



hybToGRanges = function(hybList, rna){
  seqName = rna
  hybOutput2 = hybList
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












#' Makes and ajacency matrix list (for clustering)
#'
#' Makes and ajacency matrix list (for clustering)
#' @param hybGranges list created with hybToGRanges (but just the gap section of the list)
#' @param  nucletideOrPerc measure difference by percentage or nucleotides
#' @param  cutoff The maximum difference before giving these two gaps 0
#' @return A list of Adjacancy matrices
#' @authors Jonathan Price
#' @examples 
#'  adjcancyMat =  getAdjacancyMat(hybGranges[["sample"]][["gap"]],"nucleotide", 5);
#' @export


 
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
    

    hitsWithOverlap = hits[relative_overlap <= cutoff & relative_overlap >= 0 ]
    relative_overlap = relative_overlap[relative_overlap <= cutoff  & relative_overlap >= 0 ]
  } else if(nucletideOrPerc == "perc"){

    relative_overlap = (1- (width(pintersect(x, y)) / max))

    hitsWithOverlap = hits[relative_overlap >= cutoff]
    if(length(hitsWithOverlap ) ==0 ){return(0)}else{}
    print(length(hitsWithOverlap))
    relative_overlap = relative_overlap[relative_overlap >= cutoff]

  }
  
  hitsMat = as.data.frame(hitsWithOverlap)
  hitsMat$weight = relative_overlap 


  testLong = dcast(hitsMat, queryHits ~ subjectHits, value.var = "weight")
  row.names(testLong) = testLong$queryHits
  testLong = testLong[,-1]
  testLong = as.matrix(testLong)
  testLong[is.na(testLong)] =0
  
  return(testLong)

}















#' Helper for printClusters
#'
#' 
#' @param  granges 
#' @param  indexes 
#' @param  prev 
#' @param  cluster 
#' @param  type 
#' @authors Jonathan Price
#' @examples 
#' plottingList  = printClustersFast(dir,clustering, highest_clusters, chimeraList[["sample"]][["left"]], chimeraList[["sample"]][["right"]])
#' @export

addCluster = function(granges, indexes, prev, cluster, type){
  x = granges[as.numeric(indexes)]
  x$cluster = cluster
  x$type = type
  prev= c(prev, x)
}





#' Makes a table with the coordinates of the clusters 
#'
#' 
#' @param  dir the directory that contains the *hybrids.hyb files
#' @param  clustering The output from the iGraph function cluster_walktrap for the (made with adjacency matrix input)
#' @param  highest_clusters The cluster you are interested in keeping
#' @param  left list created with hybToGRanges (but just the left section of the list)
#' @param  right list created with hybToGRanges (but just the right section of the list)
#' @return A table of clusters and coordinates
#' @authors Jonathan Price
#' @examples 
#'     plotting2 = addCluster(left, c1c1, plotting2,i  ,"left")
#' @export



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








#' Makes a table with the coordinates of the clusters 
#'
#' Does the same as printClusters but is alot fast and does not create plots
#' of each cluster 
#' @param  dir the directory that contains the *hybrids.hyb files
#' @param  clustering The output from the iGraph function cluster_walktrap for the (made with adjacency matrix input)
#' @param  highest_clusters The cluster you are interested in keeping
#' @param  left list created with hybToGRanges (but just the left section of the list)
#' @param  right list created with hybToGRanges (but just the right section of the list)
#' @return A table of clusters and coordinates
#' @authors Jonathan Price
#' @examples 
#' plottingList[[i]][[k]]  = printClustersFast(dir,clustering, highest_clusters, chimeraList[["sample"]][["left"]], chimeraList[["sample"]][["right"]])
#' @export


#helper function for printClustersFast
addCluster = function(granges, indexes, prev, cluster, type){
  x = granges[as.numeric(indexes)]
  x$cluster = cluster
  x$type = type
  prev= c(prev, x)
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









#' Plots a viewpoint of a specific region
#'
#' Takes a region of the genome and shows on a plot
#' where else that region interacts with 
#' @param  hybListSwapSars the output of swapHybs2 only
#' @param  startofChoice start of region of interest
#' @param  endofChoice The cluster you are interested in keeping
#' @param  directory where the plot should be created
#' @authors Jonathan Price
#' @examples 
#'     plotviewpoint(hybListSwapSars, 1, 100 , /path/to/plot/)
#' @export





plotviewpoint = function(hybListSwapSars,startofChoice,endofChoice, directory ){
  
 
  
  a = list()
  a[["control"]] = c(1,3)
  a[["sample"]] = c(2,4)
  
  
  df2 = data.frame(stringsAsFactors = F)
  for( i in 1:2){
    #rbind the sample
    hybListSwapSarsSample = rbind.data.frame(hybListSwapSars[[a[[i]][1]]],
                                             hybListSwapSars[[a[[i]][2]]], stringsAsFactors = F)
    
    hybListSwapSarsSampleRegion1 = hybListSwapSarsSample[((hybListSwapSarsSample$V7 < endofChoice & hybListSwapSarsSample$V7 > startofChoice) |
                                                            (hybListSwapSarsSample$V8 < endofChoice & hybListSwapSarsSample$V8 > startofChoice)|
                                                            (hybListSwapSarsSample$V8 > endofChoice & hybListSwapSarsSample$V7 < startofChoice)),]
    
    hybListSwapSarsSampleRegion2 = hybListSwapSarsSample[((hybListSwapSarsSample$V13 < endofChoice & hybListSwapSarsSample$V13 > startofChoice) |
                                                            (hybListSwapSarsSample$V14 < endofChoice & hybListSwapSarsSample$V14 > startofChoice)|
                                                            (hybListSwapSarsSample$V14 > endofChoice & hybListSwapSarsSample$V13 < startofChoice) ) , ]
    

    seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
    
    
    
    df = as.data.frame(table(unlist(seq2(from = c(hybListSwapSarsSampleRegion1$V13, hybListSwapSarsSampleRegion2$V7), to = c(hybListSwapSarsSampleRegion1$V14, hybListSwapSarsSampleRegion2$V8)))))
    
    df = df[as.numeric(as.character(df$Var1)) < startofChoice | as.numeric(as.character(df$Var1)) > endofChoice,]
    
    df$sample = names(a)[i]
    df2 = rbind.data.frame(df,df2, stringsAsFactors = F)
  }
  ggplot()+geom_line(data = df2,aes(x =as.numeric(as.character(Var1)), y = as.numeric(as.character(Freq)) ))+
    facet_grid(sample ~ . )+
    theme_classic()
  
  
  tail(df2[which(df2$Freq >300),],n = 100)
  
  pdf(paste(directory,startofChoice,"_",endofChoice,".pdf",sep = ""), height = 2.5, width = 8)
  ggplot()+geom_bar(data = df2,aes(x =as.numeric(as.character(Var1)), 
                                   y = as.numeric(as.character(Freq)) ),
                    stat = "identity")+
    facet_grid(sample~.)+
    theme_classic()
  
  dev.off()
  
  
}










#' Extracts chimeric reads for host interactions
#'
#' Extracts chimeric reads for host interactions
#' @param  TE This is the RNA of interest
#' @param  feature this is the host RNA of interest
#' @param  sampleData This is some hyb data created with readHybFiles or similar
#' @return a vector of interactions along the genome
#' @authors Jonathan Price
#' @examples 
#'          sampleVect = getVectorofInteractions(TE, feature, hybListMers[["sample"]])
#' @export






getVectorofInteractions = function(TE, feature, sampleData){
  sampleData = sampleData[sampleData$V10 == TE & sampleData$V4 == feature | sampleData$V10 == feature & sampleData$V4 == TE,]
  
  print(paste("doing ",TE,"and",feature))
  
  print("swapping")
  print(nrow(sampleData))
  if(nrow(sampleData) > 0){
    #extract the right columns
    sampleData = sampleData[,c(4,7,8,10,13,14)]
    #now swap the columns over 
    #now make sure the TE is in column 1.
    
    
    
    #Swap the columns around that do not have zika in column 2 
    sampleDatatmp1 = sampleData[sampleData$V4 == TE,]
    sampleDatatmp2 = sampleData[!(sampleData$V4 == TE),]
    sampleDatatmp2 = as.data.frame(cbind(as.character(sampleDatatmp2$V10), sampleDatatmp2$V13, sampleDatatmp2$V14, as.character(sampleDatatmp2$V4), sampleDatatmp2$V7, sampleDatatmp2$V8),stringsAsFactors = F)
    colnames(sampleDatatmp2) = c("V4","V7","V8","V10","V13","V14")
    sampleData = rbind.data.frame(sampleDatatmp1, sampleDatatmp2)
    
    #sampleData  = sampleData[sampleData$V14 > 35,]
    
    
    print("making vector")
    vect = c()
    for(i in 1:nrow(sampleData)){
      print(i)
      vect = c(vect,seq(as.numeric(sampleData$V7[i]) , as.numeric(sampleData$V8[i])))
      
    }
    return(as.data.frame(table(vect)))
  }
  
}

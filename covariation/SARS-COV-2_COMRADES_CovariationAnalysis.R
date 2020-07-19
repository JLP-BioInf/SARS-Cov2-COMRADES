# Author -  Tsveta Kamenova
# STEP-BY-STEP GUIDE TO SARS-COV-2 COMRADES CLUSTER COVARIATION ANALYSIS

library(Bios2cor)
library(seqinr)

# Step 1. Clean downloaded sequences, remove completely identical sequences and prepare files for alignment
# INPUT:
#     mySequences           a .fasta file of genome sequences
#     destinationFolder     filepath to the folder where the new file will be saved
#     newfilename           name for the new file; no extension needed
# OUTPUT:
#     a new .fasta file to be aligned 


clean.seq.set <- function(mySequences, destinationFolder, newfilename){
  
  #remove sequences containing any symbols different from A, C, T and G 
toRemove <- vector()
for(i in 1:length(mySequences)){
  for(j in 1:length(mySequences[[i]])) {
    if(mySequences[[i]][j] != "A" && mySequences[[i]][j] != "C" && mySequences[[i]][j] != "T" && 
       mySequences[[i]][j] != "G"){
      toRemove <- c(toRemove, i)
    }
  }
  print(i)
}
mySequences <- mySequences[-unique(toRemove)]

# place the reference sequence first
n <- which(names(mySequences) == "NC_045512.2")
h1 <- as.list(c2s(mySequences[[n]]))
h1 <- lapply(h1, s2c)
mySequences <- c(h1, mySequences[-n])
names(mySequences)[1] <- "SARS-CoV-2"

# remove completely identical sequences
mySequences <- lapply(mySequences, c2s)
listOfIdent <- list()
m = 1
for(i in 1:length(mySequences)){
  print(i)
  if(any(mySequences[(i+1):length(mySequences)] == mySequences[[i]])){
    g <- which(mySequences == mySequences[[i]])
    if(i == min(g)){
      listOfIdent[[m]] <- g
      m <- m+1
    } else {
      for(j in 1:length(listOfIdent)){
        if(any(listOfIdent[[j]] == i)){
          print("Added")
        }
      }
    }
  } else {
    listOfIdent[[m]] <- i
    m <- m+1
  }
  print(i)
}

vecToKeep <- sapply(listOfIdent, function(x){     #keep one sequence from each group of identical sequences
  x[1]
})
mySequences <- mySequences[vecToKeep]

filename <- paste(destinationFolder, "/", newfilename, ".fasta", sep = "")
print(filename)
write.fasta(mySequences, names = names(mySequences), file.out = filename)

}


#Step 1a (optional). Finding "seed" sequences ------------
# This code uses the kmer package to extract "seed" seqeunce identifiers
# The code uses the mbed function which computes a matrix of distances
#              between sequences in a set without alignning the sequences first.
#             Details can be found in the kmer package description, as well as in:
#                                 Blackshields, Gordon & Sievers, Fabian & Shi, Weifeng & Wilm, Andreas & Higgins, Desmond. (2010). Sequence Embedding for Fast Construction of Guide Trees for Multiple Sequence Alignment. Algorithms for molecular biology : AMB. 5. 21. 10.1186/1748-7188-5-21. 

# The following code was used for the compiling the sequences used in MSA-3-SARSrel-137seq and MSA-4-SARSrel-559seq
# For MSA-3-SARSrel-137seq, all complete non-redundant Sarbecovirus sequences (unaligned sequences from the 
#              MSA-1-SARSrel-3515seq sequence set) were clustered and seed sequences were extracted
# For  MSA-4-SARSrel-559seq, the complete non-redundant Sarbecovirus sequences were divided into 
#             seven smaller sequence sets (six 500-sequences sets and one 515-sequence set). Each set was then clustered
#             and its seed sequences were extracted and joint to form the MSA-4-SARSrel-559seq sequence set.
#
#
# INPUT:
#     mySequences           a .fasta file of sequences
#     destinationFolder     filepath to the folder where the new file will be saved
#     newfilename           name for the new file; no extension needed
# OUTPUT:
#     a new .fasta file of seed sequences to be aligned 

extract.seed.seqs <- function(mySequences, destinationFolder, newfilename){
  library(kmer)
  mySeq_mbed <- mbed(mySequences)
  seed_seq <- attr(mySeq_mbed, "seeds")
  seq_to_keep <- names(mySequences)[seed_seq]
  mySequences_new <- mySequences[seq_to_keep]
  filename <- paste(destinationFolder, "/", newfilename, ".fasta", sep = "")
  print(filename)
  write.fasta(mySequences_new, names = names(mySequences_new), file.out = filename)
}



# Step 2. MUSCLE alignment------------
#  (done externally using default MUSCLE settings and user guide instructions)
# Sequences were aligned using MUSCLE (Edgar, 2004) command line version as described in USER GUIDE at:
#             https://petrov.stanford.edu/software/src/muscle3.6/muscle3.6.html



# Step 3. Create guide file---------
#  This code creates a file of nucleotide positions in the aligned reference sequences and their corresponding
#                   positions in the non-aligned reference sequence.
# The NCBI reference SARS-CoV-2 genome (NC_045512.2) was used as reference

# INPUT:
#     alignment             a .fasta file of genome sequences
#     refseqname            sequence in the alignment to be used as reference; in this study, the
#                           NCBI reference SARS-CoV-2 genome (NC_045512.2) was used as reference
#     destinationFolder     filepath to the folder where the new file will be saved
#     newfilename           name for the new file; no extension needed
# OUTPUT:
#     a new .csv file with the nucleptide positions in the aligned reference sequences and their corresponding
#                   positions in the non-aligned reference sequence.

refseqname <- "Sars-CoV-2"

pair.checker <- function(alignment, refseqname, destinationFolder, newfilename){
  alignedSeqbyChar <- alignment[[which(names(alignment) == refseqname)]]
  
  mineToTheirs <- function(A){   # function that checks what Nt in seq correponds to position in alignment
    seqtoA <- alignedSeqbyChar[1:A]
    l1 <- length(which(seqtoA == "-"))
    intheirseqA <- as.numeric(A)-l1
    intheirseqA
  }
  
  vec1 <- 1:length(alignedSeqbyChar)
  vec2 <- alignedSeqbyChar
  vec3 <- sapply(vec1, mineToTheirs)
  
  pairChecker <- data.frame(vec1, vec2, vec3, stringsAsFactors = FALSE)
  for(i in 1:nrow(pairChecker)){
    if(pairChecker$vec2[i] == "-"){
      pairChecker$vec3[i] <- ""
    }
  }
  
  colnames(pairChecker) <- c("Position in alignment",
                             "Nucleotide", "Position in non-aligned sequence")
  
  toRem <- which(pairChecker$Nucleotide == "-")
  pair_checker <- pairChecker[ -toRem,]
  
  View(pair_checker)
  filename <- paste(destinationFolder, "/",newfilename,".csv",sep = "")
  print(filename)
  write.csv(pair_checker, file = filename, quote = F, row.names = F)
}



#Step 4. Create cluster alignments for R-scape-----------------
# This code cuts segments corresponding to the left and right side of COMRADES clusters (see text) and fuses them together for R-Scape analysis
# For each cluster side (left or right) sequences with more than 10 empty positions ("-") at the beginning
#                 of the region are discarded
#
# INPUT:
#     regionsCSV                a .csv file with coordinates of left and right sides of COMRADES clusters
#     guideCSV                  output of Step 3.: a guide .csv file with the nucleotide positions in the aligned 
#                               reference sequences and their corresponding positions in the non-aligned reference sequence.
#     myAlignment               an aligned .fasta file of sequences
#     destinationFolder         filepath to the folder where the output will be directed
#
# OUTPUT:
#     for_rscape                a folder with alignment segments corresponding of the COMRADES clusters to be analyzed with R-Scape
#     left_and_right            a folder with alignment segments corresponding of the unjoint left and right sides of COMRADES clusters
#                               (only used when there is a gap between the left and right side(gap>0))
#     removed sequences.csv     cluster numbers of the clusters which had sequences removed
#                               due to empty positions at the start of some of the extracted segnments
#     segment addresses.csv     a .csv file with paths to alignment segments generated in Step 4.


Rscape.create.regions <- function(myAlignment, regionsCSV, guideCSV, destinationFolder){
  library(Bios2cor)
  library(seqinr)
  
  for(i in 1:nrow(regionsCSV)){                         # remove spaces from cluster IDs
    regionsCSV$id[i] <- gsub(" ", "_", regionsCSV$id[i])
  }
  
  newfoldername_1 <- paste(destinationFolder, "/", "for_rscape", sep = "")
  newfoldername_2 <- paste(destinationFolder, "/", "left_and_right", sep = "")
  dir.create(newfoldername_1)
  dir.create(newfoldername_2)
  
  #   make a list every element of which corresponds to the contents of the corresponding position in the alignment
  parrot_list <- vector("list", length = length(myAlignment[[1]]))     
  for(i1 in 1:length(myAlignment[[1]])){
    for(j1 in 1:length(myAlignment)){
      parrot_list[[i1]] <- c(parrot_list[[i1]], myAlignment[[j1]][i1])
    }
    print(i1)
  }
  
  Ireduced <- vector()            # this vector will be used to store the cluster numbers of the clusters which had sequences removed
                                  #            due to empty positions at the start of some of the extracted segnments
  
  segment_addresses <- data.frame(stringsAsFactors = F)       # This dataframe will be used to store paths to the new fasta files
  alignemnt_number <- vector()
  alignemnt_address <- vector()
  
  for(i in 1:nrow(regionsCSV)){
    
    if(regionsCSV$ls[i]<regionsCSV$rs[i]){      #correct so that left side is upstream of right side
      ls <- regionsCSV$ls[i]
      le <- regionsCSV$le[i]
      rs <- regionsCSV$rs[i]
      re <- regionsCSV$re[i]
    } else {
      ls <- regionsCSV$rs[i]
      le <- regionsCSV$re[i]
      rs <- regionsCSV$ls[i]
      re <- regionsCSV$le[i]
    }
    gap <- rs - le
    print(paste("Region Set Nu.", i,  "    Left Start: ", ls, "    Left End: ", le,
                "    Right Start: ", rs, "    Right End: ", re,
                "    Gap: ", gap, sep = ""))
    
    if(gap > 0){          # if gap greater than 0 - cut out left and right side individually and join them together
      ls_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == ls)]
      le_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == le)]
      rs_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == rs)]
      re_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == re)]
      
      #cut out the LEFT side from alignment
      PL_left <- parrot_list[ls_in_alignment:le_in_alignment]
      parrot_list_rev <- vector("list", length = length(PL_left[[1]]))
      for (i1 in 1:length(PL_left[[1]])){
        for(j1 in 1:length(PL_left)){
          parrot_list_rev[[i1]] <- c(parrot_list_rev[[i1]], PL_left[[j1]][i1])
        }
      }
      
      PL_left_alignment <- parrot_list_rev
      names(PL_left_alignment) <- names(myAlignment)
      
      #remove seqs with a lot of gaps
      toRem_L <- vector("numeric")
      for(g in 1:length(PL_left_alignment)){
        if(c2s(PL_left_alignment[[g]][1:10]) == c2s(rep("-", 10))){
          toRem_L <- c(toRem_L, g)
        }
      }
      if(length(toRem_L) != 0){
        PL_left_alignment_reduced <- PL_left_alignment[-toRem_L]
        forred <- paste(as.character(i), "left")
        Ireduced <- c(Ireduced, forred)
      } else {PL_left_alignment_reduced <- PL_left_alignment}
      filename <- paste(destinationFolder, "/","left_and_right/", i, ".ID(", regionsCSV$id[i],
                      ")_REDUCED_left_", ls, "-", le, ".fasta", sep = "")
      seqinr::write.fasta(PL_left_alignment_reduced, names = names(PL_left_alignment_reduced), file.out = filename)
      
      #cut out the RIGHT side from suitable alignments
      PL_right <- parrot_list[rs_in_alignment:re_in_alignment]
      parrot_list_rev <- vector("list", length = length(PL_right[[1]]))
      for (i1 in 1:length(PL_right[[1]])){
        for(j1 in 1:length(PL_right)){
          parrot_list_rev[[i1]] <- c(parrot_list_rev[[i1]], PL_right[[j1]][i1])
        }
      }
      
      PL_right_alignment <- parrot_list_rev
      names(PL_right_alignment) <- names(myAlignment)
      
      #remove seqs with a lot of gaps
      toRem_R <- vector("numeric")
      for(g in 1:length(PL_right_alignment)){
        if(c2s(PL_right_alignment[[g]][1:10]) == c2s(rep("-", 10))){
          toRem_R <- c(toRem_R, g)
        }
      }
      length(toRem_R)
      if(length(toRem_R) != 0){
        PL_right_alignment_reduced <- PL_right_alignment[-toRem_R]
        forred <- paste(as.character(i), "right")
        Ireduced <- c(Ireduced, forred)
      } else {PL_right_alignment_reduced <- PL_right_alignment}
      filename <- paste(destinationFolder, "/", "left_and_right/", i, ".ID(", regionsCSV$id[i],
                      ")_REDUCED_right_", ls, "-", le, ".fasta", sep = "")
      seqinr::write.fasta(PL_right_alignment_reduced, names = names(PL_right_alignment_reduced), file.out = filename)
      
      
      joint_alignment <- vector("list")
      for(c in 1:length(PL_left_alignment)){
        joint_alignment[[c]] <- c(PL_left_alignment[[c]], PL_right_alignment[[which(names(PL_right_alignment) 
                                                                                    == names(PL_left_alignment)[c])]])
        names(joint_alignment)[c] <- names(PL_left_alignment)[c]
      }
     
      # join the reduced alignments
      reduced_joint_alignment <- vector("list")
      m = 1
      for(c in 1:length(PL_left_alignment_reduced)){
        if(any(names(PL_right_alignment_reduced) == names(PL_left_alignment_reduced)[c])){
          reduced_joint_alignment[[m]] <- c(PL_left_alignment_reduced[[c]], 
                                            PL_right_alignment_reduced[[which(names(PL_right_alignment_reduced) 
                                                                              == names(PL_left_alignment_reduced)[c])]])
          names(reduced_joint_alignment)[m] <- names(PL_left_alignment_reduced)[c]
          m = m+1
        }
      }
     
      filename <- paste(destinationFolder, "/", "for_rscape/", i, ".ID(", regionsCSV$id[i],
                      ")_REDUCED_(", ls, "-", le, ")_(", rs, "-", re, ").fasta", sep = "")
      seqinr::write.fasta(reduced_joint_alignment, names = names(reduced_joint_alignment), file.out = filename)
      alignemnt_address[i] <- filename
      alignemnt_number[i] <- i
      
    } else {         # if gap smaller than 0 (sides overlap) take the whole region between the start of the left side and the end of the right side
      ls_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == ls)]
      re_in_alignment <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == re)]
      
      #cut region out from suitable alignments
      PL_LandR <- parrot_list[ls_in_alignment:re_in_alignment]
      
      parrot_list_rev <- vector("list", length = length(PL_LandR[[1]]))
      for (i1 in 1:length(PL_LandR[[1]])){
        for(j1 in 1:length(PL_LandR)){
          parrot_list_rev[[i1]] <- c(parrot_list_rev[[i1]], PL_LandR[[j1]][i1])
        }
      }
      
      PL_alignment <- parrot_list_rev
      names(PL_alignment) <- names(myAlignment)
      
      # removed seqs with gaps
      toRem_LandR<- vector("numeric")
      for(g in 1:length(PL_alignment)){
        if(c2s(PL_alignment[[g]][1:10]) == c2s(rep("-", 10))){
          toRem_LandR <- c(toRem_LandR, g)
        }
      }
      length(toRem_LandR)
      if(length(toRem_LandR) != 0){
        PL_alignment_reduced <- PL_alignment[-toRem_LandR]
        forred <- as.character(i)
        Ireduced <- c(Ireduced, forred)
      } else {PL_alignment_reduced <- PL_alignment}
    
      filename <- paste(destinationFolder, "/", "for_rscape/", i, ".ID(", regionsCSV$id[i],
                      ")_REDUCED_(", ls, "-", le, ")_(", rs, "-", re, ").fasta", sep = "")
      seqinr::write.fasta(PL_alignment_reduced, names = names(PL_alignment_reduced), filename)
      alignemnt_address[i] <- filename
      alignemnt_number[i] <- i
    }
  }

  filename <- paste(destinationFolder, "/", "removed sequences.csv", sep = "")        # Save csv of cluster numbers of clusters with removed sequences
  write.csv(Ireduced, file = filename, row.names = FALSE, quote = FALSE)
  
  segment_addresses <- cbind(alignemnt_number, alignemnt_address)
  colnames(segment_addresses) <- c("Nu.", "Path")
  filename <- paste(destinationFolder, "/", "segment addresses.csv", sep = "")
  write.csv(segment_addresses, file = filename, quote = FALSE, row.names = FALSE)
  
}




# Step 5. R-scape analysis-------------
# The alignemnt segments in the for_rscape folder generated in Step 4. were analysed with R-Scape (Rivas et al., 2017)
#           as a command line tool. The results were compiled in a .txt file which will be analyzed in Step 6.
# The pipeline used for the R-Scape analysis of the alignment segments
#           in the for_rscape folder is given in the CommandLine_RScape_pipe.txt file (included).

# Step 6. Analysis of R-scape results---------------
# This function analyzes the .txt file of R-Scape results.
# INPUT:
#     myResults                 .txt file of R-Scape results for all alignment segments generated in Step 4.
#     regionsCSV                a .csv file with coordinates of left and right sides of COMRADES clusters used 
#                               in for generation of the alignment segments
#     guideCSV                  output of Step 3.: a guide .csv file with the nucleotide positions in the aligned 
#                               reference sequences and their corresponding positions in the non-aligned reference sequence.
#                               Also used in Step 4. for generation of the alignment segments.
#                                
#     segmentAddresses          the segment addresses.csv generated in Step 4.
#     destinationFolder         filepath to the folder where the output files will be directed
#     newfilename               common identifier name for the output files
#
# OUTPUT:
#     newfilename - Results from all clusters.csv                        a .csv file with the output from each clusters
#     newfilename - Table of suggested co-varying positions.csv          a .csv file with all position combinations suggested to covary by R-Scape    
#                                                                        (This file also includes position that cannot form base pairs)
#     newfilename - BPs >90% of seqs.csv                                 a .csv file with unique potential base pairs only


RScape.results.analysis <- function(myResults, regionsCSV, guideCSV, segmentAddresses, destinationFolder, newfilename){
  colnames(myResults) <- "res"
  
  rows <- which(myResults$res == myResults$res[4])
  rows_odd <- rows[seq(1, length(rows), 2)]
  rows_even <-  rows[seq(2, length(rows), 2)]
  MSA_names <- as.character(myResults$res[sapply(rows_odd, function(x){
    y <- x+1
    y
  })])
  MSA_names <- sapply(MSA_names, as.character)
  #View(MSA_names)
  MSA_names_short <- sapply(MSA_names, function(x){
    y <- unlist(strsplit(x, " "))[3]
    y
  })
  
  result_rows <- sapply(rows_even, function(x){
    y <- x+1
    y
  })
  rscape_output <- myResults$res[result_rows]
  rscape_output <- sapply(rscape_output, as.character)
  
  # make a table with the output of R-scape from each alignment
  ALL_RES_TABLE <- data.frame(MSA_names_short, rscape_output, stringsAsFactors = FALSE)
  colnames(ALL_RES_TABLE) <- c("MSA name", "Output")
  for(i in 1:nrow(ALL_RES_TABLE)){
    if(as.character(ALL_RES_TABLE$Output[i]) == "~"){
      ALL_RES_TABLE$Output[i] <- as.character("------------BASE-PAIRS IDENTIFIED - SEE TABLE OF SUGGESTED CO-VARYING POSITIONS------------")
    }
  }
  filename <- paste(destinationFolder, "/", newfilename, " - Results from all clusters.csv", sep = "")
  print(filename)
  write.csv(ALL_RES_TABLE, file = filename, quote = F, row.names = F)
  
  
  # Compile full results table
  pos_with_sug_BP <- which(ALL_RES_TABLE$Output == "------------BASE-PAIRS IDENTIFIED - SEE TABLE OF SUGGESTED CO-VARYING POSITIONS------------")
  length(pos_with_sug_BP)
  
  tableBPsIdentified <- data.frame(stringsAsFactors = FALSE)
  ClusterNum <- vector("numeric")
  clusterID <- vector("character")
  LS <- vector()
  LE <- vector()
  RS <- vector()
  RE <- vector()
  gap <- vector()
  pos1_in_alignment <- vector()
  pos2_in_alignment <- vector()
  score <- vector()
  e_Value <- vector()
  subst <- vector()
  powerstat <- vector()
  
  j = 0
  m = 1
  for(i in 1:length(pos_with_sug_BP)){
    workingA <- ALL_RES_TABLE$`MSA name`[pos_with_sug_BP[i]]
    a <- as.numeric(unlist(strsplit(workingA, "ID"))[1])
    ClusterNum[m] <- a
    clusterID[m] <- as.character(regionsCSV$id[a])
    if(regionsCSV$ls[a]<regionsCSV$rs[a]){
      LS[m] <- as.numeric(regionsCSV$ls[a])
      LE[m] <- as.numeric(regionsCSV$le[a])
      RS[m] <- as.numeric(regionsCSV$rs[a])
      RE[m] <- as.numeric(regionsCSV$re[a])
    } else {
      LS[m] <- as.numeric(regionsCSV$rs[a])
      LE[m] <- as.numeric(regionsCSV$re[a])
      RS[m] <- as.numeric(regionsCSV$ls[a])
      RE[m] <- as.numeric(regionsCSV$le[a])
    }
    gap[m] <- RS[m] - LE[m]
    b <- as.numeric(result_rows[pos_with_sug_BP[i]])
    c <- as.character(myResults$res[seq(b, b+8, 1)])
    pos1_in_alignment[m] <- as.numeric(gsub(" ", "", c[2]))
    pos2_in_alignment[m] <- as.numeric(gsub(" ", "", c[3]))
    score[m] <- as.numeric(c[4]) 
    e_Value[m] <- as.numeric(c[5])
    subst[m] <- as.numeric(c[6])
    powerstat[m] <- as.numeric(c[8])
    if(c[9] == "~"){
      j = 1
    } else {j = 0}
    while(j == 1){
      m <- m +1
      ClusterNum[m] <- a
      clusterID[m] <- as.character(regionsCSV$id[a])
      if(regionsCSV$ls[a]<regionsCSV$rs[a]){
        LS[m] <- as.numeric(regionsCSV$ls[a])
        LE[m] <- as.numeric(regionsCSV$le[a])
        RS[m] <- as.numeric(regionsCSV$rs[a])
        RE[m] <- as.numeric(regionsCSV$re[a])
      } else {
        LS[m] <- as.numeric(regionsCSV$rs[a])
        LE[m] <- as.numeric(regionsCSV$re[a])
        RS[m] <- as.numeric(regionsCSV$ls[a])
        RE[m] <- as.numeric(regionsCSV$le[a])
      }
      gap[m] <- RS[m] - LE[m]
      b <- b+8
      c <- as.character(myResults$res[seq(b, b+8, 1)])
      pos1_in_alignment[m] <- as.numeric(gsub(" ", "", c[2]))
      pos2_in_alignment[m] <- as.numeric(gsub(" ", "", c[3]))
      score[m] <- as.numeric(c[4]) 
      e_Value[m] <- as.numeric(c[5])
      subst[m] <- as.numeric(c[6])
      powerstat[m] <- as.numeric(c[8])
      if(c[9] == "~"){
        j = 1
      } else {
        j = 0
      }
    }
    m = m+1
  }
  tableBPsIdentified <- data.frame(ClusterNum, clusterID, LS, LE, RS, RE, gap, pos1_in_alignment, pos2_in_alignment, score, e_Value, subst, powerstat, stringsAsFactors = FALSE)
  colnames(tableBPsIdentified) <- c("Cluster number" ,"Cluster ID", "LS", "LE", "RS", "RE","Gap", "Pos 1 in alignemnt", "Pos 2 in alignment", "Score", "E-value", "Substitutions", "Power")
  
  pos1_in_seq <- vector()
  pos1_in_seq_nt <- vector()
  pos2_in_seq_nt <- vector()
  pos2_in_seq <- vector()
  for(i in 1:nrow(tableBPsIdentified)){
    pos1_in_alignment <- tableBPsIdentified$`Pos 1 in alignemnt`[i]
    pos2_in_alignment <- tableBPsIdentified$`Pos 2 in alignment`[i]
    lsAli <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == tableBPsIdentified$LS[i])]
    leAli <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == tableBPsIdentified$LE[i])]
    rsAli <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == tableBPsIdentified$RS[i])]
    reAli <- guideCSV$Position.in.alignment[which(guideCSV$Position.in.non.aligned.sequence == tableBPsIdentified$RE[i])]
    gap <- tableBPsIdentified$Gap[i]
    if(gap >= 0){
      pos_indeces <- c(seq(lsAli,leAli,1), seq(rsAli,reAli,1))
    } else {
      pos_indeces <- seq(lsAli,reAli,1)
    }
    
    pos1_index <- pos_indeces[tableBPsIdentified$`Pos 1 in alignemnt`[i]]
    pos2_index <- pos_indeces[tableBPsIdentified$`Pos 2 in alignment`[i]]
    
    if(any(guideCSV$Position.in.alignment == pos1_index)){
      row1 <- which(guideCSV$Position.in.alignment == pos1_index)
      pos1_in_seq[i] <- guideCSV$Position.in.non.aligned.sequence[row1]
      pos1_in_seq_nt[i] <- as.character(guideCSV$Nucleotide[row1])
    }
    if(any(guideCSV$Position.in.alignment == pos2_index)){
      row2 <- which(guideCSV$Position.in.alignment == pos2_index)
      pos2_in_seq[i] <- guideCSV$Position.in.non.aligned.sequence[row2]
      pos2_in_seq_nt[i] <- as.character(guideCSV$Nucleotide[row2])
    }
  }
  
  tableBPsIdentified <- cbind(tableBPsIdentified, pos1_in_seq, pos1_in_seq_nt, pos2_in_seq, pos2_in_seq_nt)
  colnames(tableBPsIdentified)[14:17] <- c("Position 1 in non-aligned Sequence", "Nucleotide 1 in non-aligned Sequence", 
                                           "Position 2 in non-aligned Sequence", "Nucleotide 2 in non-aligned Sequence")
  
  # Add columns of percentage of base-pairing nucleotides at each positions combination, as well as a column of frequencies of the various nucleotide 
  # combinations observed at the respective positions combination
  suggested_pos <- tableBPsIdentified
  bp_percent <- vector("numeric")
  content_perc <- vector("character")
  for(i in 1:nrow(suggested_pos)){
    print(paste("Calculating nucleotide combination frequencies - position pair No.", i, sep = ""))
    clust_num <- suggested_pos$`Cluster number`[i]
    alignment <- import.fasta(segmentAddresses[clust_num,2])
    parrot_list <- vector("list", length = length(alignment[[1]]))
    for (i1 in 1:length(alignment[[1]])){
      for(j1 in 1:length(alignment)){
        parrot_list[[i1]] <- c(parrot_list[[i1]], alignment[[j1]][i1])
      }
    }
    
    pos1_in_segment <- suggested_pos$`Pos 1 in alignemnt`[i]
    pos2_in_segment <- suggested_pos$`Pos 2 in alignment`[i]
    
    pos1_full <- parrot_list[pos1_in_segment]
    pos2_full <- parrot_list[pos2_in_segment]
    
    pangolin <- data.frame(pos1_full, pos2_full, stringsAsFactors = F)
    colnames(pangolin) <- c("Nuc1", "Nuc2")
    
    nt_pairs <- vector("character", length = nrow(pangolin))
    for(j in 1:nrow(pangolin)){
      nt_pairs[j] <- paste(pangolin$Nuc1[j], pangolin$Nuc2[j], sep = "--")
    }
    pangolin <- cbind(pangolin, nt_pairs)
    colnames(pangolin)[3] <- "Nuc Pair"
    
    bear <- sort(summary(factor(pangolin$`Nuc Pair`)), decreasing = TRUE)
    q <- vector()
    for(k in 1:length(bear)){
      if(k != length(bear)){
        q_name <- names(bear)[k]
        q_perc <- (bear[k] / (sum(bear))) * 100
        q_perc <- round(q_perc, 2)
        q[k] <- paste(as.character(q_name), " (", as.character(q_perc), "); ", sep = "")
      } else {
        q_name <- names(bear)[k]
        q_perc <- (bear[k] / (sum(bear))) * 100
        q_perc <- round(q_perc, 2)
        q[k] <- paste(q_name, " (",q_perc, ")", sep = "")
      }
    }
    
    bp <- 0
    non_bp <- 0
    for(k in 1:length(pangolin$`Nuc Pair`)){
      x <- pangolin$`Nuc Pair`[k]
      if(x == "A--T" | x == "T--A" | x == "G--T" | x == "T--G" | x == "G--C" | x == "C--G"){
        bp <- bp+1
      } else { non_bp <- non_bp+1 }
    }
    bp_percent[i] <- (bp/(bp+non_bp))*100
    bp_percent[i] <- round(bp_percent[i], 2)
    
    bear_united <- as.character(q[1])
    for(m in 2:length(q)){
      bear_united <- paste(bear_united, as.character(q[m]), sep = )
    }
    content_perc[i] <- bear_united 
  }
  suggested_pos <- cbind(suggested_pos, bp_percent, content_perc)
  a <- ncol(suggested_pos)
  colnames(suggested_pos)[(a-1):a] <- c("Per.cent.base.pairing.nt.combination", "Position.nt.pair.composition")
  
  filename <- paste(destinationFolder, "/", newfilename, " - Table of suggested co-varying positions.csv",  sep = "")
  print(filename)
  write.csv(suggested_pos, file =  filename, quote = F, row.names = F)
  
  
  # Remove positions combinations that CANNOT form a base pair in more than 105 of the sequences in the alignment
  reduced_sugg_pos <- subset(suggested_pos, suggested_pos$Per.cent.base.pairing.nt.combination >= 90)
  
  lastcol <- vector(length = nrow(reduced_sugg_pos))
  for(o in 1:length(lastcol)){
    lastcol[o] <- reduced_sugg_pos$`Position 2 in non-aligned Sequence`[o] - reduced_sugg_pos$`Position 1 in non-aligned Sequence`[o]
  }
  reduced_sugg_pos <- cbind(reduced_sugg_pos, lastcol)
  reduced_sugg_pos <- subset(reduced_sugg_pos, reduced_sugg_pos$lastcol > 3)
  reduced_sugg_pos <- reduced_sugg_pos[-20]
  nrow(reduced_sugg_pos)
  mean(reduced_sugg_pos$`E-value`)
  
  uniquePairs <- vector()
  toKeep_2 <- vector()
  for(v in 1:nrow(reduced_sugg_pos)){
    pos1 <- reduced_sugg_pos$`Position 1 in non-aligned Sequence`[v]
    pos2 <- reduced_sugg_pos$`Position 2 in non-aligned Sequence`[v]
    pos_pair <- paste(as.character(pos1), "--", as.character(pos2), sep = "")
    if(any(uniquePairs == pos_pair)){
    } else {
      uniquePairs <- c(uniquePairs, pos_pair)
      toKeep_2 <- c(toKeep_2, v)
    }
  }
  reduced_sugg_pos <- reduced_sugg_pos[toKeep_2,]
  uniquePairs
  nrow(reduced_sugg_pos)
  View(reduced_sugg_pos)
  filename <- paste(destinationFolder, "/", newfilename, " - BPs >90% of seqs.csv",  sep = "")
  print(filename)
  write.csv(reduced_sugg_pos, file =  filename, quote = F, row.names = F)
}


# Step 7. Create a summary table with results from all alignments and samples-------------
# The following code combines all potential base pairs from the R-Scape analysis of four alignments (MSA-1-SARSrel-3515seq, MSA-2-SARSrel-3515seq, MSA-3-SARSrel-137seq and MSA-4-SARSrel-559seq, see text)
#                       and two SARS-CoV-2 COMRADES sample (Sample 1 and Sample 2)
# The full lists of suggested position combination are available in the AlignmentAndSample_full_results_files folder
# In this form, the code reproduces Supplementary Table 2 and saves the output file in the desired folder (destionationFolder)
# 
# INPUT:
#     destinationFolder                                path to the folder where the output file will be directed
#
# OUTPUT:
#     SARS-CoV-2_COMRADES_Covariation_SUMMARY.csv      a list of unique position combinations identified in any of the four alignemnts and two samples analyzed


rm(list = ls())



combine.result.files <- function(destinationFolder){

  #Results files:
  al137_S1 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS137_clustS1 - Table of suggested co-varying positions.csv")
  al137_S2 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS137_clustS2 - Table of suggested co-varying positions.csv")
  al559_S1 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS559_clustS1 - Table of suggested co-varying positions.csv")
  al559_S2 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS559_clustS2 - Table of suggested co-varying positions.csv")
  al3515_S1 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS3515_clustS1 - Table of suggested co-varying positions.csv")
  al3515_S2 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS3515_clustS2 - Table of suggested co-varying positions.csv")
  al3515_S1_maxit2 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS3515_muscleMaxIter2_clustS1 - Table of suggested co-varying positions.csv")
  al3515_S2_maxit2 <- read.csv("./3.AlignmentAndSample_full_results_files/SARS3515_muscleMaxIter2_clustS2 - Table of suggested co-varying positions.csv")
  
  
  combine.pairs <- function(x){                # this function adds two extra columns ("Position (pair)" and "Nt [nucleotide] pair") that facilitate subsequent analysis
  pos_pair <- vector()
  nt_pair <- vector()
  tableBPsIdentified <- x
  for(i in 1:nrow(tableBPsIdentified)){
    pos1 <- tableBPsIdentified$Position.1.in.non.aligned.Sequence[i]
    pos2 <- tableBPsIdentified$Position.2.in.non.aligned.Sequence[i]
    pos_pair[i] <- paste(as.character(pos1), "--", as.character(pos2), sep = "")
    
    nt1 <- tableBPsIdentified$Nucleotide.1.in.non.aligned.Sequence[i]
    nt2 <- tableBPsIdentified$Nucleotide.2.in.non.aligned.Sequence[i]
    nt_pair[i] <- paste(nt1, "--", nt2, sep = "")
  }
  tableBPsIdentified <- cbind(tableBPsIdentified, pos_pair, nt_pair)
  colnames(tableBPsIdentified)[(ncol(tableBPsIdentified)-1): (ncol(tableBPsIdentified))] <- c("Positions (pair)", ("Nt pair"))
  x <- tableBPsIdentified
  x
}


al137_S1 <- combine.pairs(al137_S1)             #add the two extra columns to each alignment
al137_S2 <- combine.pairs(al137_S2)
al559_S1 <- combine.pairs(al559_S1)
al559_S2 <- combine.pairs(al559_S2)
al3515_S1 <- combine.pairs(al3515_S1)
al3515_S2 <- combine.pairs(al3515_S2)
al3515_S1_maxit2 <- combine.pairs(al3515_S1_maxit2)
al3515_S2_maxit2 <- combine.pairs(al3515_S2_maxit2)


al137_S1 <- cbind(al137_S1, "al137_S1")             #add an identifier to each column
al137_S2 <- cbind(al137_S2, "al137_S2")
al559_S1 <- cbind(al559_S1, "al559_S1")
al559_S2 <- cbind(al559_S2, "al559_S2")
al3515_S1 <- cbind(al3515_S1, "al3515_S1")
al3515_S2 <- cbind(al3515_S2, "al3515_S2")
al3515_S1_maxit2 <- cbind(al3515_S1_maxit2, "al3515_S1_maxit2")
al3515_S2_maxit2 <- cbind(al3515_S2_maxit2, "al3515_S2_maxit2")


colnames(al3515_S1)[22] <- "identifieralignemnt"

colnames(al137_S1) <- colnames(al3515_S1)
colnames(al137_S2) <- colnames(al3515_S1)
colnames(al559_S1) <- colnames(al3515_S1)
colnames(al559_S2) <- colnames(al3515_S1)
colnames(al3515_S1) <- colnames(al3515_S1)
colnames(al3515_S2) <- colnames(al3515_S1)
colnames(al3515_S1_maxit2) <- colnames(al3515_S1)
colnames(al3515_S2_maxit2) <- colnames(al3515_S1)

#create the joint table

summarytable <- data.frame(stringsAsFactors = FALSE)
summarytable <- rbind(summarytable, al137_S1, al137_S2, al559_S1, al559_S2,
                      al3515_S1, al3515_S2, al3515_S1_maxit2, al3515_S2_maxit2)
summarytable <- subset(summarytable, summarytable$Per.cent.base.pairing.nt.combination > 90)

lastcol <- vector(length = nrow(summarytable))
for(j in 1:length(lastcol)){
  lastcol[j] <- summarytable$Position.2.in.non.aligned.Sequence[j] - summarytable$Position.1.in.non.aligned.Sequence[j]
}
summarytable <- cbind(summarytable, lastcol)
summarytable <- subset(summarytable, summarytable$lastcol > 3)
ncol(summarytable)
summarytable <- summarytable[ ,-23]


uniquePairs <- vector()
toKeep_2 <- vector()
for(i in 1:nrow(summarytable)){
  pos1 <- summarytable$Position.1.in.non.aligned.Sequence[i]
  pos2 <- summarytable$Position.2.in.non.aligned.Sequence[i]
  pos_pair <- paste(as.character(pos1), "--", as.character(pos2), sep = "")
  if(any(uniquePairs == pos_pair)){
  } else {
    uniquePairs <- c(uniquePairs, pos_pair)
    toKeep_2 <- c(toKeep_2, i)
  }
}
summarytable_reduced <- summarytable[toKeep_2,]

for(i in 1:nrow(summarytable_reduced)){
  myBP <- summarytable_reduced$`Positions (pair)`[i]
  working_table <- subset(summarytable, summarytable$`Positions (pair)` == myBP)
  mean_E_val <- mean(working_table$E.value)
  summarytable_reduced$E.value[i] <- mean_E_val
}
which(colnames(summarytable_reduced) == "E.value")
colnames(summarytable_reduced)[11] <- "Average.E.value"


#View(summarytable_reduced)

Occ.SARS137_S1 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS137_S2 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS559_S1 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS559_S2 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS3515_S1 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS3515_S2 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS3515_S1_maxiter2 <- vector("character", length = nrow(summarytable_reduced))
Occ.SARS3515_S2_maxiter2 <- vector("character", length = nrow(summarytable_reduced))

summarytablewithoccurences <- cbind(summarytable_reduced, Occ.SARS3515_S1, Occ.SARS3515_S2,
                                    Occ.SARS3515_S1_maxiter2, Occ.SARS3515_S2_maxiter2,
                                    Occ.SARS137_S1, Occ.SARS137_S2, 
                                    Occ.SARS559_S1, Occ.SARS559_S2)


# find where each bp is found
for(i in 1:nrow(summarytablewithoccurences)){
  print(i)
  working_table <- subset(summarytable, summarytable$`Positions (pair)` == summarytablewithoccurences$`Positions (pair)`[i])
  m <- unique(working_table$identifieralignemnt)
  
  #al137_S1
  if(any(m == "al137_S1")){
    summarytablewithoccurences$Occ.SARS137_S1[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS137_S1[i] <- ""
  }
  
  #al137_S2
  if(any(m == "al137_S2")){
    summarytablewithoccurences$Occ.SARS137_S2[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS137_S2[i] <- ""
  }
  
  #al559_S1
  if(any(m == "al559_S1")){
    summarytablewithoccurences$Occ.SARS559_S1[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS559_S1[i] <- ""
  }
  
  #al559_S2
  if(any(m == "al559_S2")){
    summarytablewithoccurences$Occ.SARS559_S2[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS559_S2[i] <- ""
  }
  
  #al3515_S1
  if(any(m == "al3515_S1")){
    summarytablewithoccurences$Occ.SARS3515_S1[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS3515_S1[i] <- ""
  }
  
  #al3515_S2
  if(any(m == "al3515_S2")){
    summarytablewithoccurences$Occ.SARS3515_S2[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS3515_S2[i] <- ""
  }
  
  #al3515_S1_maxit2
  if(any(m == "al3515_S1_maxit2")){
    summarytablewithoccurences$Occ.SARS3515_S1_maxiter2[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS3515_S1_maxiter2[i] <- ""
  }
  
  #al3515_S2_maxit2
  if(any(m == "al3515_S2_maxit2")){
    summarytablewithoccurences$Occ.SARS3515_S2_maxiter2[i] <- "YES"
  } else {
    summarytablewithoccurences$Occ.SARS3515_S2_maxiter2[i] <- ""
  }
}

which(colnames(summarytablewithoccurences) == "identifieralignemnt")
length(which(summarytablewithoccurences$Occ.SARS559_S2 == "YES"))
summarytablewithoccurences <- summarytablewithoccurences[ , -22]
s <- summarytablewithoccurences
clean_summarytablewithoccurences <- cbind(s$Position.1.in.non.aligned.Sequence, s$Nucleotide.1.in.non.aligned.Sequence, 
                                           s$Position.2.in.non.aligned.Sequence, s$Nucleotide.2.in.non.aligned.Sequence,
                                           s$Average.E.value, s$Occ.SARS3515_S1, s$Occ.SARS3515_S2, s$Occ.SARS3515_S1_maxiter2, s$Occ.SARS3515_S2_maxiter2,
                                           s$Occ.SARS137_S1, s$Occ.SARS137_S2, s$Occ.SARS559_S1, s$Occ.SARS559_S2)
colnames(clean_summarytablewithoccurences) <- c("Position.1.in.non.aligned.Sequence", "Nucleotide.1.in.non.aligned.Sequence", 
                                                "Position.2.in.non.aligned.Sequence", "Nucleotide.2.in.non.aligned.Sequence", colnames(s)[11],
                                                      "MSA-1-SARSrel-3515seq(Sample1)", "MSA-1-SARSrel-3515seq(Sample2)",
                                                      "MSA-2-SARSrel-3515seq(Sample1)", "MSA-2-SARSrel-3515seq(Sample2)",
                                                      "MSA-3-SARSrel-137seq(Sample1)", "MSA-3-SARSrel-137seq(Sample2)",
                                                      "MSA-4-SARSrel-559seq(Sample1)", "MSA-4-SARSrel-559seq(Sample2)")
filename <- paste(destinationFolder, "/SARS-CoV-2_COMRADES_Covariation_SUMMARY.csv",  sep = "")
print(filename)
write.csv(clean_summarytablewithoccurences, file = filename,
          quote = FALSE, row.names = FALSE)
}




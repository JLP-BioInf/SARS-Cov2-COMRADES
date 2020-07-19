################################################################
# This code is specific to the project:                        #
# The short- and long-range RNA-RNA Interactome of SARS-CoV-2  #
# It creates contact heatmaps from raw data                    #
################################################################




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

source("COMRADESfunctions.R")

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
# Make a matrix for the contact maps
# and plot (inlcuding the line traces )
##############################

library(RColorBrewer)
hybListMersSwap = swapHybs2(hybListMers, IDMers)
hybListSarsSwap = swapHybs(hybListSars, IDSars)

matrixListMers = getMatrices(hybListMersSwap, IDMers)
matrixListSars = getMatrices(hybListSarsSwap, IDSars)



plotMatrixList(matrixListMers,  "./figures/contactMaps/mers", TE, hybDirMers)
plotMatrixList(matrixListSars,  "./figures/contactMaps/sars", TE, hybDirSars)


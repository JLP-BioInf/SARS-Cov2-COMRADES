################################################################
# This code is specific to the project:                        #
# The short- and long-range RNA-RNA Interactome of SARS-CoV-2  #
#  It plot specific viewpoints for SARS and MERS               #
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
# get the Hyb files swapped
##############################
hybListSwapMers = swapHybs2(hybListMers, IDMers)
hybListSarsSwap = swapHybs2(hybListSars, IDSars)

##############################
# choose a region to plot the viepoint of
##############################


# 5 prime UTR viewpoint
plotviewpoint(hybListSarsSwap, 1, 400, ".")
plotviewpoint(hybListSwapMers, 1, 400, ".")

# 3 prime UTR viewpoint
plotviewpoint(hybListSarsSwap, 29500, 30000, ".")




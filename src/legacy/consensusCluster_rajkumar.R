################################################################################### METADATA ####
#
# consensusCluster_rajkumar.R
# Version: 1.0.1 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Perform unsupervised clustering on imputed counts.
#
###################################################################################### LEGAL ####
#
# Copyright © 2020 King's College London
#
# This work is licensed under the Creative Commons Attribution 4.0 International Licence. To view
# a copy of this license,  visit http://creativecommons.org/licences/by/4.0/  or send a letter to
# Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
#
# Permission is hereby granted,  free of charge, to any person  obtaining a copy of this software
# and  associated  documentation  files  (the  "Software"),  to  deal  in  the  Software  without
# restriction,  including without  limitation the  rights to use,  copy, modify,  merge, publish,
# distribute, and/or sell copies  of the Software, and to permit persons  to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above  copyright notice  and this  permission notice  shall be included  with all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED  "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR IMPLIED, INCLUDING
# BUT NOT  LIMITED TO THE WARRANTIES  OF MERCHANTABILITY,  FITNESS FOR A  PARTICULAR PURPOSE, AND
# NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS  OR COPYRIGHT HOLDERS BE  LIABLE FOR ANY CLAIM,
# DAMAGES, OR  OTHER LIABILITY,  WHETHER IN AN  ACTION OF CONTRACT,  TORT, OR  OTHERWISE, ARISING
# FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF OR DEALING IN THE SOFTWARE.
#
########################################################################## FILE REQUIREMENTS ####
#
# [1] combineCounts/Data/metadata.xlsx
#     Curated metadata file.
#
##################################################################################### INPUTS ####
#
# [1] combineCounts/Results/combineCounts.tsv
#     Combined normalised counts files in tsv format.
#
#################################################################################### OUTPUTS ####
#
# [1] console figure viewer
#     ConsensusClusterPlus plots.
#
# [2] supplementary_rajkumar/Results/Clustering/
#     Cluster membership details.
#
############################################################################# INITIALISATION ####
require(ConsensusClusterPlus)
require(data.table)
require(dplyr)
require(foreach)
require(Rmagic)
require(readxl)
require(here)

## File input

infile <- "combineCounts/Results/combineCounts.tsv"                        # Counts matrix
inExcel <- "combineCounts/Data/metadata.xlsx"                                                # colData
outdir <- "supplementary_rajkumar/Results/Clustering/"

dir.create(outdir,
           showWarnings = F)

## Read file

colData <- read_excel(inExcel, sheet = 1)                                                                   # Read colData
data <- fread(infile)                                                                                       # Read counts
genes <- data[[1]]                                                                                          # Record genes in counts file

## Select samples

colData_select <- colData %>%                                                                               # Select which samples you want to accept
    filter(grepl("rajkumar", Source)) %>%
    filter(Phenotype == "Parkinson's disease" | Phenotype == "Alzheimer's disease") %>%
    dplyr::select(Colname) %>%
    unlist() %>%
    as.character()

data_select <- data %>%
    dplyr::select(all_of(colData_select))

data_select <- data_select %>% as.matrix()
rownames(data_select) <- genes

## Prepare data matrix for MAGIC

data_select_MAGIC <- t(data_select)                                                                         # Transpose (now genes are on columns, samples are on rows)
genes_all_na <- which(is.na(colMeans(data_select_MAGIC, na.rm = T))) %>% unlist() %>% as.numeric()          # Find the columns (genes) which are all NA
if(length(genes_all_na) != 0) {
    data_select_MAGIC <- data_select_MAGIC[,-genes_all_na]                                                  # Remove the columns with all NA if there are any
}
data_select_MAGIC[is.na(data_select_MAGIC)] <- 0                                                            # Set all remaining NAs to zero (MAGIC will impute zero values)
genes_all_zero <- which(colSums(data_select_MAGIC) == 0) %>% unlist() %>% as.numeric()                      # Find the columns which are all zero
if(length(genes_all_zero) != 0) {
    data_select_MAGIC <- data_select_MAGIC[,-genes_all_zero]                                                # Remove the columns with all zero (as these cannot be imputed) if there are any
}

## Perform MAGIC impute

data_select_MAGIC <- magic(data = data_select_MAGIC, npca = 50)
data_select_MAGIC_CC <- data_select_MAGIC$result %>% as.matrix()
data_select_MAGIC_CC <- t(data_select_MAGIC_CC)

## Normalise

mads <- apply(X = data_select_MAGIC_CC, MARGIN = 1, FUN = mad)                                              # Normalise
data_select_MAGIC_CC <- data_select_MAGIC_CC[rev(x = order(mads))[1:5000],]
data_select_MAGIC_CC <- sweep(x = data_select_MAGIC_CC, MARGIN = 1,
                              STATS = apply(X = data_select_MAGIC_CC, MARGIN = 1, FUN = median, na.rm=T))

## Consensus cluster

wd <- here()
setwd(outdir)
res <- ConsensusClusterPlus(d = data_select_MAGIC_CC, maxK = 7, reps = 1000, distance = "pearson", clusterAlg = "hc" , plot = "png")
setwd(wd)

## Record consensus classes

nClusters <- 3                                                                                              # Set number of clusters to accept

classes <- foreach (i = 1:nClusters,                                                                        # Record members of each cluster
                    .final = function(i) setNames(i, 1:nClusters)) %do% {
                        which(res[[nClusters]]$consensusClass == i) %>% names()
                    }

for (i in 1:length(classes)) {
    fwrite(x = as.list(classes[[i]]),
           file = paste0(outdir, "cluster_", i, ".tsv"),
           sep = "\t")
}


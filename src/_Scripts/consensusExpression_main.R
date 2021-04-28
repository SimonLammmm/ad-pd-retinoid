################################################################################### METADATA ####
#
# consensusExpression_main.R
# Version: 1.0.2 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Determine mean expression counts within clusters.
#
###################################################################################### LEGAL ####
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
# [2] combineCounts/Results/combineCounts_entrez.tsv
#     Combined normalised counts files in tsv format.
#
# [3] main/Results/Clustering/metadataWithClusters.tsv
#     Metadata with added cluster membership column in tsv format.
#
#################################################################################### OUTPUTS ####
#
# [1] main/Results/ConsensusExpression/cluster_X_consensus_expression.tsv
#     Mean expression values for cluster X
#
# [2] main/Results/ConsensusExpression/countsConsensusEntrez.tsv
#     Mean expression values per cluster (Entrez format)
#
############################################################################# INITIALISATION ####
require(data.table)
require(dplyr)
require(foreach)

inClusters <- "main/Results/Clustering/metadataWithClusters.tsv"
inCounts <- "combineCounts/Results/combineCounts.tsv"
inCountsEntrez <- "combineCounts/Results/combineCounts_entrez.tsv"
outdir <- "main/Results/ConsensusExpression/"

dir.create(outdir, showWarnings = F)

## Determine cluster members and number of clusters

clusters <- fread(inClusters)

nClusters <- max(clusters$Cluster)

clusterMembers <- foreach(i = 1:nClusters,
                          .final = function(i) setNames(i, 1:nClusters)) %do% {
                              clusters %>% filter(Cluster == i) %>% dplyr::select(Colname) %>% unlist() %>% as.character()
                          }

## Read counts file

countsFull <- fread(inCounts)
countsFullEntrez <- fread(inCountsEntrez)

## Subset counts file for cluster members

getCountsByCluster <- function(clusterMembers, countsFull) {

    countsByCluster <- foreach(i = clusterMembers,
                               .final = function(i) setNames(i, names(clusterMembers))) %do% {
                                   countsFull %>% dplyr::select(i)
                               }

    return(countsByCluster)

}

countsByCluster <- getCountsByCluster(clusterMembers, countsFull)

## Determine mean expression values per cluster and write to file (one file per cluster)

getCountsConsensus <- function(countsByCluster, countsFull) {

    countsConsensus <- foreach(i = countsByCluster,
                               .final = function(i) setNames(i, names(countsByCluster))) %do% {
                                   cbind(countsFull[,1],
                                         rowMeans(i, na.rm = T))
                               }

    return(countsConsensus)

}

countsConsensus <- getCountsConsensus(countsByCluster, countsFull)

##

for (i in 1:length(countsConsensus)) {
    fwrite(countsConsensus[[i]], paste0(outdir, "cluster_", i, "_consensus_expression.tsv"),
           sep = "\t")
}


## Determine mean expression values per cluster and write to file (one file total, entrez genes) - for CellFie

countsConsensusEntrez <- getCountsByCluster(clusterMembers, countsFullEntrez) %>% getCountsConsensus(countsFullEntrez)

countsConsensusEntrezUnlist <- countsConsensusEntrez[[1]][,1]

for (i in 1:length(countsConsensusEntrez)) {
    countsConsensusEntrezUnlist <- cbind(countsConsensusEntrezUnlist, countsConsensusEntrez[[i]][,2])
}

names(countsConsensusEntrezUnlist) <- c("genes", 1:(length(countsConsensusEntrezUnlist)-1))

##

fwrite(countsConsensusEntrezUnlist, paste0(outdir, "countsConsensusEntrez.tsv"),
       sep = "\t")

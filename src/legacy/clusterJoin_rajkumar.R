################################################################################### METADATA ####
#
# clusterJoin_rajkumar.R
# Version: 1.0.1 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Add a column to the metadata denoting cluster membership from unsupervised clustering.
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
# [1] supplementary_rajkumar/Results/Clustering/
#     Cluster membership details.
#
#################################################################################### OUTPUTS ####
#
# [1] supplementary_rajkumar/Results/Clustering/metadataWithClusters.tsv
#     Metadata with added cluster membership column in tsv format.
#
############################################################################# INITIALISATION ####
require(dplyr)
require(data.table)
require(foreach)
require(readxl)

infile_colData <- "combineCounts/Data/metadata.xlsx"
indir_clusters <- "supplementary_rajkumar/Results/Clustering/"

outfile <- "supplementary_rajkumar/Results/Clustering/metadataWithClusters.tsv"

infile_clusters <- list.files(indir_clusters, full.names = T, pattern = "^cluster_\\d+.*\\.tsv")
infile_clusters_names <- list.files(indir_clusters, pattern = "^cluster_\\d+.*\\.tsv")

colData <- read_xlsx(infile_colData, sheet = 1)

clustersList <- foreach (i = 1:length(infile_clusters),
                         .final = function(i) setNames(i, infile_clusters_names)) %do% {
                             tibble(Colname = fread(infile_clusters[i]) %>% names(),
                                    Cluster = i) %>%
                                 left_join(colData, by = "Colname")
                         }

clustersOut <- tibble()
for (clust in clustersList) {
    clustersOut <- rbind(clustersOut, clust)
}

fwrite(clustersOut,
       file = outfile,
       sep = "\t")

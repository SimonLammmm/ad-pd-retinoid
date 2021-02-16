################################################################################### METADATA ####
#
# Vis_rajkumar.R
# Version: 1.0.3 (2020-10-16)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Exploratory data visualisations.
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
##################################################################################### INPUTS ####
#
# [1] vis/Data/combineCounts_impute.tsv
#     Combined normalised imputed counts in tsv format.
#
# [2] combineCounts/Data/metadata.xlsx
#     Curated metadata file.
#
# [3] supplementary_rajkumar/Results/Clustering/metadataWithClusters.tsv
#     Metadata with added cluster membership column in tsv format.
#
#################################################################################### OUTPUTS ####
#
# [1] supplementary_rajkumar/Results/vis/
#     Visualisations.
#
############################################################################# INITIALISATION ####
require(dplyr)
require(ggplot2)
require(Rtsne)
require(data.table)
require(readxl)
#require(ggbiplot)
require(umap)

## IO

incounts <- "vis/Data/combineCounts_impute.tsv"
cts3_MAGIC_vis <- as.matrix(fread(incounts))

incolData <- "combineCounts/Data/metadata.xlsx"
colData <- read_excel(incolData, sheet = 1)

inclusters <- "supplementary_rajkumar/Results/Clustering/metadataWithClusters.tsv"
clusters <- fread(inclusters)

outplot <- "supplementary_rajkumar/Results/vis/"

set.seed(31415)                                                                                             # For a reproducible plot

####################################################################################### MAIN ####

## t-SNE

cts3_tsne <- Rtsne(cts3_MAGIC_vis, check_duplicates = F)
cts3_tsne2 <- data.frame(tSNE1 = cts3_tsne$Y[,1],
                         tSNE2 = cts3_tsne$Y[,2]) %>%
    cbind(colData) %>%
    full_join(clusters %>% dplyr::select(Colname, Cluster), by = "Colname")

png(paste0(outplot, "tSNE_by_phenotype.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = Phenotype)) + geom_point() +  ggtitle("tSNE by phenotype")
dev.off()

png(paste0(outplot, "tSNE_by_source.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = Source)) + geom_point() + ggtitle("tSNE by source")
dev.off()

png(paste0(outplot, "tSNE_by_cluster.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = factor(Cluster))) + geom_point() + ggtitle("tSNE by cluster")
dev.off()

png(paste0(outplot, "tSNE_by_organ.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = Region1)) + geom_point() + ggtitle("tSNE by organ")
dev.off()

png(paste0(outplot, "tSNE_by_brain_region_1.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = Region2)) + geom_point() + ggtitle("tSNE by brain region")
dev.off()

png(paste0(outplot, "tSNE_by_brain_region_2.png"), width = 800, height = 600)
ggplot(cts3_tsne2, aes(x = tSNE1, y = tSNE2, colour = Region3)) + geom_point() + ggtitle("tSNE by brain region")
dev.off()


## PCA

cts3_pca <- prcomp(cts3_MAGIC_vis)
cts3_pca2 <- data.frame(pca1 = cts3_pca$x[,1],
                        pca2 = cts3_pca$x[,2]) %>%
    cbind(colData) %>%
    full_join(clusters %>% dplyr::select(Colname, Cluster), by = "Colname")

png(paste0(outplot, "pca_by_phenotype.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Phenotype)) + geom_point() +  ggtitle("PCA by phenotype")
dev.off()

png(paste0(outplot, "pca_by_source.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Source)) + geom_point() + ggtitle("PCA by source")
dev.off()

png(paste0(outplot, "pca_by_source_zoom.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Source)) + geom_point() + ggtitle("PCA by source") + xlim(-75,70) + ylim(-7.5,5)
dev.off()

png(paste0(outplot, "pca_by_cluster.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = factor(Cluster))) + geom_point() + ggtitle("PCA by cluster")
dev.off()

png(paste0(outplot, "pca_by_organ.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Region1)) + geom_point() + ggtitle("PCA by organ")
dev.off()

png(paste0(outplot, "pca_by_brain_region_1.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Region2)) + geom_point() + ggtitle("PCA by brain region")
dev.off()

png(paste0(outplot, "pca_by_brain_region_2.png"), width = 800, height = 600)
ggplot(cts3_pca2, aes(x = pca1, y = pca2, colour = Region3)) + geom_point() + ggtitle("PCA by brain region")
dev.off()


## UMAP

cts3_umap <- umap(cts3_MAGIC_vis)
cts3_umap2 <- data.frame(umap1 = cts3_umap$layout[,1],
                         umap2 = cts3_umap$layout[,2]) %>%
    cbind(colData) %>%
    full_join(clusters %>% dplyr::select(Colname, Cluster), by = "Colname")

png(paste0(outplot, "umap_by_phenotype.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Phenotype)) + geom_point() +  ggtitle("UMAP by phenotype")
dev.off()

png(paste0(outplot, "umap_by_source.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Source)) + geom_point() + ggtitle("UMAP by source")
dev.off()

png(paste0(outplot, "umap_by_source_zoom.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Source)) + geom_point() + ggtitle("UMAP by source") + xlim(-15, 25) + ylim(-30, 30)
dev.off()

png(paste0(outplot, "umap_by_cluster.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = factor(Cluster))) + geom_point() + ggtitle("UMAP by cluster")
dev.off()

png(paste0(outplot, "umap_by_organ.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Region1)) + geom_point() + ggtitle("UMAP by organ")
dev.off()

png(paste0(outplot, "umap_by_brain_region_1.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Region2)) + geom_point() + ggtitle("UMAP by brain region")
dev.off()

png(paste0(outplot, "umap_by_brain_region_2.png"), width = 800, height = 600)
ggplot(cts3_umap2, aes(x = umap1, y = umap2, colour = Region3)) + geom_point() + ggtitle("UMAP by brain region")
dev.off()

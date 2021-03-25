################################################################################### METADATA ####
#
# DESeq2_main.R
# Version: 1.0.5 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Perform differential gene expression analysis with DESeq2.
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
# [1] main/Results/Clustering/
#     Cluster membership details.
#
# [2] combineCounts/Results/combineCounts.tsv
#     Combined normalised counts file.
#
#################################################################################### OUTPUTS ####
#
# [1] main/Results/DESeq2/
#     DESeq2 results.
#
############################################################################# INITIALISATION ####
require(data.table)
require(DESeq2)
require(dplyr)
require(foreach)
require(ggplot2)
require(readxl)
require(scales)

indir <- "main/Results/Clustering/"                                                       # cluster membership file dir
indir_filefilter <- "cluster_"                                                                              # filter keyword for cluster membership files
inExcel <- "combineCounts/Data/metadata.xlsx"                                                # colData
inCounts <- "combineCounts/Results/combineCounts.tsv"                      # Counts matrix
user_alpha <- 1e-10                                                                                         # Set alpha for determining significant DEGs
user_l2fc_thresh <- 0                                                                                       # Set l2fc threshold for determining significant DEGs
outdir <- "main/Results/DESeq2/"                                                                   # Output dir for saved files

dir.create(outdir,
           showWarnings = F)

################################################################################## FUNCTIONS ####

see <- function(x) {                                                                                        # Helper function to tabulate DESeq results
    tibble(gene = x@rownames,
           baseMean = x$baseMean,
           l2fc = x$log2FoldChange,
           l2fc_se = x$lfcSE,
           stat = x$stat,
           pval = x$pvalue,
           padj = x$padj)
}

reverselog_trans <- function(base = exp(1)) {                                                               # Inverse log scale for ggplot2
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

####################################################################################### MAIN ####

# File I/O

infiles <- list.files(path = indir,                                                                         # List clusters files
                      pattern = indir_filefilter)

infiles_full <- list.files(path = indir,
                           pattern = indir_filefilter,
                           full.names = T)

data <- foreach (file = infiles_full,                                                                       # List clusters and memberships
                 .final = function(file) setNames(file, infiles)) %do% {
                     fread(file) %>% names()
                 }

data_unlist <- data %>% unlist() %>% as.character()

colData <- read_excel(inExcel, sheet = 1)                                                                   # Read colData

counts <- fread(inCounts)                                                                                   # Read counts matrix
genes <- counts[[1]]                                                                                        # Record gene names
counts <- counts[,-1]                                                                                       # Remove gene column from actual matrix
rownames(counts) <- genes                                                                                   # Add back gene names as row headers

# Data selection

colData_cluster <- colData %>%                                                                              # Add cluster data to colData
    filter(Colname %in% data_unlist) %>%
    filter(!(Colname %in% names(counts)[colMeans(counts, na.rm = T) %>% unlist() %>% is.nan()])) %>%        # Drop columns with all NA counts
    mutate(Cluster = case_when(Colname %in% data[[1]] ~ 1,                                                  ### ADD MORE CASES IF YOU HAVE MORE CLUSTERS
                               Colname %in% data[[2]] ~ 2,                                                  ### ADD MORE CASES IF YOU HAVE MORE CLUSTERS
                               Colname %in% data[[3]] ~ 3,                                                  ### ADD MORE CASES IF YOU HAVE MORE CLUSTERS
                               Colname %in% data[[4]] ~ 4))

colData_cluster$Cluster <- as.factor(colData_cluster$Cluster)

colData_select <- colData_cluster %>%                                                                       # Filter colData for only relevant samples (those in any cluster)
    dplyr::select(Colname) %>%
    unlist() %>%
    as.character()

counts_select <- counts %>%                                                                                 # Filter counts for only relevant samples
    dplyr::select(colData_select) %>%
    as.matrix() %>%
    exp() %>%
    sapply(as.integer)

# Data formatting

excess_NA_threshold <- dim(counts_select)[2] * 0.4                                                          # Set the threshold for removing genes if too many NAs are present
genes_idx_excess_NA <- which(rowCounts(x = counts_select, value = NA) >= excess_NA_threshold)               # Remove genes with too many NA
if (length(genes_idx_excess_NA) != 0) {
    counts_select <- counts_select[-genes_idx_excess_NA,]
    genes <- genes[-genes_idx_excess_NA]
}

#genes_idx_all_NA <- which(is.na(rowMeans(counts_select, na.rm = T)))                                       # Remove genes with all NA
#if (length(genes_idx_all_NA) != 0) {
#    counts_select <- counts_select[-genes_idx_all_NA,]
#    genes <- genes[-genes_idx_all_NA]
#}

genes_idx_all_zero <- which(rowMeans(counts_select, na.rm = T) == 0) %>% unlist() %>% as.numeric()          # Remove genes with all zero
if (length(genes_idx_all_zero) != 0) {
    counts_select <- counts_select[-genes_idx_all_zero,]
    genes <- genes[-genes_idx_all_zero]
}

#counts_select <- sapply(X = counts_select, FUN = as.integer)                                               # Convert counts matrix to integer
#counts_select <- exp(counts_select)                                                                        # Adjust counts so that no values are negative
counts_select[is.na(counts_select)] <- 0                                                                    # NA to zero (DESeq assumes not detected)

rownames(counts_select) <- genes                                                                            # Restore gene names

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts_select,
                              colData = colData_cluster,
                              design = ~ Cluster)

sizeFactors(dds) <- 1                                                                                       # Tell DESeq that the counts are already normalised

dds <- DESeq(dds)
# cluster b (differential), cluster a (reference)
a <- 4

for (b in 1:(a-1)) {
    
    res <- results(dds,
                   contrast = c("Cluster", b, a))
    
    
    resData <- see(res)
    
    # Volcano plot
    
    resData <- resData %>%
        mutate(colour = case_when(abs(l2fc) < user_l2fc_thresh | padj > user_alpha ~ "grey",
                                  abs(l2fc) >= user_l2fc_thresh & padj <= user_alpha ~ "red")) %>%
        arrange(padj)
    
    ggplot(data = resData,
           mapping = aes(x = l2fc, y = padj, colour = colour, alpha = colour)) +
        geom_point() +
        scale_alpha_manual(values = c(0.1, 0.3, 0)) +
        scale_colour_manual(values = c("#AAAAAA", "#FF0000", "#000000")) +
        scale_y_continuous(trans = reverselog_trans(10)) +
        #geom_vline(xintercept = -user_l2fc_thresh,
        #           col = "black",
        #           linetype = "dotted",
        #           size = 0.7) +
        #geom_vline(xintercept = user_l2fc_thresh,
        #           col = "black",
        #           linetype = "dotted",
        #           size = 0.7) +
        geom_hline(yintercept = user_alpha,
                   col = "black",
                   linetype = "dotted",
                   size = 0.7) +
        theme(legend.position = "none") +
        labs(title = paste0("Volcano plot: cluster ", b, " versus cluster ", a),
             x = "Log 2 fold change",
             y = "Adjusted p value")
    ggsave(paste0(outdir, "DESeq2_c", b, "_vs_c", a, ".png"), device = "png")
    
    
    # Significant DEGs
    
    resData_sig <- resData %>%
        filter(colour == "red")
    
    # Save to file
    
    fwrite(x = resData,
           file = paste0(outdir,
                         "allDEG_c",
                         b,
                         "_vs_c",
                         a,
                         ".tsv"),
           sep = "\t")
    
    fwrite(x = resData_sig,
           file = paste0(outdir,
                         "sigDEG_c",
                         b,
                         "_vs_c",
                         a,
                         ".tsv"),
           sep = "\t")
    
}

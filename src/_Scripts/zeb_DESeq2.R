################################################################################### METADATA ####
#
# zeb_DESeq2.R
# Version: 1.0.2 (2020-04-28)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Programmatically run and analyse multiple pairwise differential expression analysis with DESeq2.
# In zebrafish expression analysis pipeline.
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
# [1] zeb/Data/finalCts/finalCts.se.Rdata
#     Cleaned and combined counts and adjacent metadata ready for DESeq2 analysis.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Results/DESeq2
#     DESeq2 results.
#
################################################################################### OVERVIEW ####
#
# 0. Initialise environment
#    a. Attach packages
#    b. Register parallelisation
#    c. Define user variables
#    d. Define functions
# 1. Load summarizedExperiments from file
# 2. Do DESeq
# 3. Gather results
# 4. Plot DESeq results
#    a. MA plot
#    b. Volcano plot
# 5. Identify significant genes
# 6. Clean up environment
#
############################################################################# INITIALISATION ####

### Initialisation of packages
require(SummarizedExperiment)
require(DESeq2)
require(doParallel)
require(foreach)
require(dplyr)
require(biomaRt)
require(ggplot2)
require(scales)
require(gridExtra)
require(data.table)

### Initialisation of parallel environment
registerDoParallel(detectCores() - 2)

userSE <- "zeb/Data/finalCts/finalCts.se.Rdata"
user_alpha <- 0.1                                                       # Alpha for significance testing
#user_l2fc_thresh <- 2                                                   # Threshold for significant l2fc
user_out <- "zeb/Results/DESeq2/"

################################################################################## FUNCTIONS ####

cf <- function(x, y, z) {                                               # Helper function to get model variables and levels
    sort(as.character(unique(x[[y]][[z]])))
}

see <- function(x) {                                                    # Helper function to tabulate DESeq results
    tibble(gene = x@rownames,
           baseMean = x$baseMean,
           l2fc = x$log2FoldChange,
           l2fc_se = x$lfcSE,
           stat = x$stat,
           pval = x$pvalue,
           padj = x$padj)
}

getBioMart <- function(genes) {                                         # Convert gene id to symbol with biomaRt
    getBM(attributes = c("ensembl_transcript_id",
                         "hgnc_symbol",
                         "description"),
          filters = "ensembl_transcript_id",
          values = genes,
          mart = useMart("ensembl",
                         host = "https://www.ensembl.org",
                         dataset = "drerio_gene_ensembl"))              # Querying zebrafish genes
}

reverselog_trans <- function(base = exp(1)) {                           # Inverse log scale for ggplot2
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

masterPlot <- function(dataset) {                                       # Master ggplot2 volcano plotter
    foreach (plot = names(dataset)) %dopar% {
        plotName <- strsplit(x = plot, split = "\\.")
        plotdata <- see(dataset[[plot]]) %>%
            mutate(colour = case_when(abs(l2fc) < 2 | padj > user_alpha ~ "grey",
                                      abs(l2fc) >= 2 & padj <= user_alpha ~ "red"))
        ggplot(data = plotdata,
               mapping = aes(x = l2fc, y = padj, colour = colour, alpha = colour)) +
            geom_point() +
            scale_alpha_manual(values = c(0.1, 0.3, 0)) +
            scale_colour_manual(values = c("#AAAAAA", "#FF0000", "#000000")) +
            scale_y_continuous(trans = reverselog_trans(10)) +
            #geom_vline(xintercept = -user_l2fc_thresh,
            #           col = "black",
            #           linetype = "dotted",
            #           size = 0.7) +
#            geom_vline(xintercept = user_l2fc_thresh,
#                       col = "black",
#                       linetype = "dotted",
#                       size = 0.7) +
            geom_hline(yintercept = user_alpha,
                       col = "black",
                       linetype = "dotted",
                       size = 0.7) +
            theme(legend.position = "none") +
            labs(title = paste(#"Volcano plot: ",
                               plotName[[1]][1],
                               ", ",
                               plotName[[1]][2],
                               " vs. ",
                               plotName[[1]][3],
                               sep = ""),
                 xlab = "Log 2 fold change",
                 ylab = "Adjusted p value")
    }
}


####################################################################################### MAIN ####

### 1. Load summarizedExperiments from file
load(file = userSE)                                                     # Import summarizedExperiments

dds <- list()
for (exp in names(se)) {                                                # Import as DESeq dataset
    print(exp)
    try({
        dds[[exp]] <- DESeqDataSet(se = se[[exp]],                      # Try to import with design ~ part
                                   design = ~ part)
    },
    silent = T)
    try({
        dds[[exp]] <- DESeqDataSet(se = se[[exp]],                      # Try to import with design ~ tert1
                                   design = ~ tert1)
    },
    silent = T)
    try({
        dds[[exp]] <- DESeqDataSet(se = se[[exp]],                      # Try to import with design ~ part + tert1
                                   design = ~ part + tert1 + part:tert1)
    },
    silent = T)
}

### 2. Do DESeq
deseq <- foreach (exp = names(dds),
                  .final = function(exp) setNames(exp, names(dds))) %dopar% {    # Do DESeq on each se
                      DESeq(dds[[exp]])
                  }

### 3. Gather results
resultsPart <- foreach (exp = names(deseq),                             # Get DESeq results by part
                        .final = function(exp) setNames(exp, names(deseq))) %dopar% {
                            foreach (exppart1 = cf(deseq, exp, "part"),
                                     .final = function(exppart1) setNames(exppart1, cf(deseq, exp, "part"))) %dopar% {
                                         foreach (exppart2 = cf(deseq, exp, "part"),
                                                  .final = function(exppart2) setNames(exppart2, cf(deseq, exp, "part"))) %dopar% {
                                                      if (exppart1 < exppart2) {
                                                          results(deseq[[exp]],
                                                                  contrast = c("part", exppart2, exppart1))
                                                      }
                                                  }
                                     }
                        }
resultsGen <- foreach (exp = names(deseq),                              # Get DESeq results by genotype
                       .final = function(exp) setNames(exp, names(deseq))) %dopar% {
                           foreach (expgen1 = cf(deseq, exp, "tert1"),
                                    .final = function(expgen1) setNames(expgen1, cf(deseq, exp, "tert1"))) %dopar% {
                                        foreach (expgen2 = cf(deseq, exp, "tert1"),
                                                 .final = function(expgen2) setNames(expgen2, cf(deseq, exp, "tert1"))) %dopar% {
                                                     if (expgen1 < expgen2) {
                                                         results(deseq[[exp]],
                                                                 contrast = c("tert1", expgen2, expgen1))
                                                     }
                                                 }
                                    }
                       }

resultsPart <- unlist(resultsPart)                                      # Collapse lists
resultsGen <- unlist(resultsGen)
parts <- cf(deseq, "Full", "part")
gens <- cf(deseq, "Full", "tert1")

### 4. Plot DESeq results
## 4a. MA plot
#par(mfrow = c(5, 3))                                                    # Plot DESeq by genotype
#par(mar = c(1,1,1,1))
#for (row in c("Full", parts)) {
#    for (vs1 in gens) {
#        for (vs2 in gens) {
#            if (vs1 < vs2) {
#                h <- paste(row, vs1, vs2, sep = ".")
#                j <- paste(vs1, "vs", vs2)
#                plotMA(resultsGen[[h]], ylim = c(-10,10))
#                title(main = paste0(row, ": ", j))
#            }
#        }
#    }
#}
#
#par(mfrow = c(4, 6))
#par(mar = c(1,1,1,1))
#for (row in c("Full", gens)) {                                          # Plot DESeq by part
#    for (vs1 in parts) {
#        for (vs2 in parts) {
#            if (vs1 < vs2) {
#                h <- paste(row, vs1, vs2, sep = ".")
#                j <- paste(vs1, "vs", vs2)
#                plotMA(resultsPart[[h]], ylim = c(-10,10))
#                title(main = paste0(row, ": ", j))
#            }
#        }
#    }
#}

## 4b. Volcano plot
dir.create(paste0(user_out, "Plots"), showWarnings = F)

plotsByPart <- masterPlot(resultsPart)                                  # Create plots
png(file = paste0(user_out, "Plots/Volcano_plot_comparing_genotypes.png"),
    height = 1500,
    width = 2500)
print({
	do.call(grid.arrange, c(plotsByPart, nrow = 4))                         # Plot on a grid
})
dev.off()

plotsByGen <- masterPlot(resultsGen)
png(file = paste0(user_out, "Plots/Volcano_plot_comparing_parts.png"),
    height = 1500,
    width = 1000)
print({
	do.call(grid.arrange, c(plotsByGen, nrow = 5))
})
dev.off()

### 5. Identify significant genes
allGenesByGenotype <- foreach (res = names(resultsGen),                 # Get all genes by genotype
                               .final = function(res) setNames(res, names(resultsGen))) %dopar% {
                                   see(resultsGen[[res]]) %>%
                                       arrange(padj) %>%
                                       mutate(subset = strsplit(x = res, split = "\\.")[[1]][1],
                                              differential = strsplit(x = res, split = "\\.")[[1]][2],
                                              reference = strsplit(x = res, split = "\\.")[[1]][3])
                               }
allGenesByGenotype <- do.call(rbind, allGenesByGenotype)                # Collapse into a single tibble

sigGenesByGenotype <- allGenesByGenotype %>%                            # Get significant genes by genotype
    filter(padj <= user_alpha)# %>%
  #  filter(abs(l2fc) >= user_l2fc_thresh)
biomaRtGenesByGenotype <- getBioMart(genes = sigGenesByGenotype$gene)   # Get gene annotation from biomaRt for significant genes only
sigGenesByGenotype <- sigGenesByGenotype %>%
    full_join(biomaRtGenesByGenotype,
              by = c("gene" = "ensembl_transcript_id")) %>%
    unique()

allGenesByPart <- foreach (res = names(resultsPart),                    # Get all genes by part
                           .final = function(res) setNames(res, names(resultsPart))) %dopar% {
                               see(resultsPart[[res]]) %>%
                                   filter(padj < user_alpha) %>%
                                  # filter(abs(l2fc) >= user_l2fc_thresh) %>%
                                   arrange(padj) %>%
                                   mutate(subset = strsplit(x = res, split = "\\.")[[1]][1],
                                          differential = strsplit(x = res, split = "\\.")[[1]][2],
                                          reference = strsplit(x = res, split = "\\.")[[1]][3])
                           }
allGenesByPart <- do.call(rbind, allGenesByPart)                        # Collapse into a single tibble

sigGenesByPart <- allGenesByPart %>%                                    # Get significant genes by genotype
    filter(padj <= user_alpha) #%>%
    #filter(abs(l2fc) >= user_l2fc_thresh)
biomaRtGenesByPart <- getBioMart(genes = sigGenesByPart$gene)   # Get gene annotation from biomaRt for significant genes only
sigGenesByPart <- sigGenesByPart %>%
    full_join(biomaRtGenesByPart,
              by = c("gene" = "ensembl_transcript_id")) %>%
    unique()

allSummaryByGenotype <- allGenesByGenotype %>%                            # Count number of total genes across all parts
    group_by(subset, differential, reference) %>% summarise(genes = n())

allSummaryByPart <- allGenesByPart %>%                                    # and across all genotypes
    group_by(subset, differential, reference) %>% summarise(genes = n())

sigSummaryByGenotype <- sigGenesByGenotype %>%                            # Count number of significant genes across all parts
    group_by(subset, differential, reference) %>% summarise(genes = n())

sigSummaryByPart <- sigGenesByPart %>%                                    # and across all genotypes
    group_by(subset, differential, reference) %>% summarise(genes = n())

dir.create(paste0(user_out, "sig"),
           showWarnings = F)
dir.create(paste0(user_out, "all"),
           showWarnings = F)

fwrite(x = allGenesByGenotype,                                         # Write to file
       file = paste0(user_out, "all/zeballGenes_ByGenotype.tsv"))
fwrite(x = allGenesByPart,
       file = paste0(user_out, "all/zeballGenes_ByPart.tsv"))
fwrite(x = allSummaryByGenotype,
       file = paste0(user_out, "all/zeballGenes_ByGenotype_summary.tsv"))
fwrite(x = allSummaryByPart,
       file = paste0(user_out, "all/zeballGenes_ByPart_summary.tsv"))
fwrite(x = sigGenesByGenotype,
       file = paste0(user_out, "all/zebsigGenes_ByGenotype.tsv"))
fwrite(x = sigGenesByPart,
       file = paste0(user_out, "all/zebsigGenes_ByPart.tsv"))
fwrite(x = sigSummaryByGenotype,
       file = paste0(user_out, "all/zebsigGenes_ByGenotype_summary.tsv"))
fwrite(x = sigSummaryByPart,
       file = paste0(user_out, "all/zebsigGenes_ByPart_summary.tsv"))

### 6. Clean up environment
stopImplicitCluster()

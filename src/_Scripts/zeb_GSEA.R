################################################################################### METADATA ####
#
# zeb_GSEA.R
# Version: 1.0.0 (2020-02-26)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Programmatically  perform piano  gene set  enrichment analysis  on a set of  DESeq2 results. In
# zebrafish expression analysis pipeline.
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
# [1] zeb/Results/DESeq2
#     DESeq2 results.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Results/GSEA
#     GSEA results.
#
################################################################################### OVERVIEW ####
#
# 0. Initialise environment
#    a. Attach packages
#    b. Register parallelisation
#    c. Define functions
# 1. Wrangle input
#    a. Read input file
#    b. Subset input file by experiment
# 2. Generate inputs for piano
#    a. Determine l2fc directions
#    b. Record p values (not p adj)
#    c. Query BioMart for gene sets
# 3. Run piano GSEA
# 4. Record results
# 5. Clean up environment
#
############################################################################# INITIALISATION ####

### Initialisation of packages
require(data.table)
require(dplyr)
require(tidyr)
require(doParallel)
require(foreach)
require(biomaRt)
require(piano)

### Initialisation of parallel environment
registerDoParallel(detectCores() - 2)

### User variables
userGenes <- "zeb/Results/DESeq2/all/zeballGenes_ByGenotype.tsv"        # Results file from DESeq2 pipeline
userSets <- "zeb/Data/GSEA/biomaRtForGSC.tsv"                                 # Gene sets for piano GSEA
userOut <- "zeb/Results/GSEA/"

################################################################################## FUNCTIONS ####

getBioMart <- function(genes) {                                         # Convert gene id to symbol with biomaRt
    getBM(attributes = c("ensembl_transcript_id",
                         "name_1006"),
          filters = "ensembl_transcript_id",
          values = genes,
          mart = useMart("ensembl",
                         dataset = "drerio_gene_ensembl"))              # Querying zebrafish genes

}

### Get biomaRt results
fetchBM <- function(x,                                                               # Character or vector to submit to biomaRt
                    file = NULL                                                      # File to check for reading or writing
) {
    # If file is specified, see if biomaRt results already exist, and if so, load them
    if (!is.null(file)) {                                                            # Check if file is specified
        try(
            suppressWarnings(res <- read.delim(file = file,                          # Try to load the file
                                               colClasses = "character")),
            silent = TRUE
        )
    }
    # If biomaRt results don't exist locally, fetch them from biomaRt and save them to file if specified
    while (!exists("res")) {                                                         # Check if biomaRt results were loaded successfully
        try(                                                                         # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("ensembl_transcript_id",
                                        "name_1006"),
                         filters = "ensembl_transcript_id",
                         values = x,
                         mart = useMart("ensembl",
                                        dataset = "drerio_gene_ensembl")))           # Querying zebrafish genes
        if (exists("res") & (!is.null(file))) {                                      # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

####################################################################################### MAIN ####

## 1. Wrangle input
# 1a. Read input file
genesRaw <- fread(userGenes)                                            # Read in DEGs from file
contrasts <- genesRaw %>%
    group_by(subset, differential, reference) %>%
    summarise() %>%
    mutate(name = paste(subset, differential, reference, sep = "_"))

# 1b. Subset input file by experiment
genes <- foreach (i = 1:nrow(contrasts),                                # Subset DEGs by contrast
                  .final = function(i) setNames(i, contrasts$name)) %do% {
    genesRaw %>%
        filter(subset==as.character(contrasts[i,1])) %>%
        filter(differential==as.character(contrasts[i,2])) %>%
        filter(reference==as.character(contrasts[i,3])) %>%
        filter(!is.na(padj)) %>%
        mutate(dir = sign(l2fc)) %>%
        dplyr::select(gene, pval, dir)
                  }

## 2. Generate inputs for piano
# 2a. Determine l2fc directions
dir = list()                                                            # Collect directions for each DEG
for (n in names(genes)) {
    dir[[n]] <- genes[[n]]$dir
    names(dir[[n]]) <- genes[[n]]$gene
}

# 2b. Record p values (not p adj)
pval = list()                                                           # Collect p values for each DEG
for (n in names(genes)) {
    pval[[n]] <- genes[[n]]$pval
    names(pval[[n]]) <- genes[[n]]$gene
}

# 2c. Query BioMart for gene sets
sets <- fetchBM(unique(genesRaw$gene), userSets) %>%                    # Generate gene sets
    filter(name_1006 != "") %>%
    loadGSC()

## 3. Run piano GSEA
gsaRes <- foreach (n = names(genes),
                   .final = function(n) setNames(n, names(genes))) %do% {
                       runGSA(geneLevelStats = pval[[n]],
                              directions = dir[[n]],
                              gsc = sets)
                   }

## 4. Record results
# 4a. Tables

foreach (n = names(genes)) %do% {
    fwrite(x = GSAsummaryTable(gsaRes[[n]]),
           file = paste(userOut,
                        n,
                        "_GSEAsummary.tsv",
                        sep = ""))
    fwrite(x = data.frame(names = names(gsaRes[[n]]$gsc), genes = as.character(gsaRes[[n]]$gsc)),
           file = paste(userOut,
                        n,
                        "_gsc_genes.tsv",
                        sep = ""))
}
# 4b. Heatmaps
heights <- c(1400, 1800, 2100, 2800, 3700,
             3100, 20200, 1200, 2600, 2700,
             2700, 1800, 2800, 2200, 1500)
foreach (i = 1:length(names(genes))) %do% {
    png(file = paste(userOut,
                     names(genes)[i],
                     "_GSEAsummary.png",
                     sep = ""),
        width = 940,
        height = heights[i])
    GSAheatmap(gsaRes[[i]],
               cutoff = 20,
               adjusted = T,
               ncharLabel = 64,
               cellnote = "none",
               colorkey = F,
               cex = 2)
    dev.off()
}
#png(file = paste0(userOut, "GSEAlegend.png"),
#    width = 900,
#    height = 900)
#GSAheatmap(gsaRes[[1]],
#           cutoff = 1,
#           adjusted = T,
#           ncharLabel = 64,
#           cellnote = "none",
#           cex = 2)
#dev.off()

## 5. Clean up environment
stopImplicitCluster()

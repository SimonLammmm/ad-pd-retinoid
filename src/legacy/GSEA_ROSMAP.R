################################################################################### METADATA ####
#
# GSEA_ROSMAP.R
# Version: 1.0.6 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Perform gene set enrichment analysis with piano.
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
# [1] supplementary_ROSMAP/Results/DESeq2/
#     DESeq2 results.
#
# [-] supplementary_ROSMAP/Data/GSEA/BioMartforGSC.tsv
#     KEGG dictionary for building  the gene set collections.  Optional file, will be dynamically
#     created if not supplied.
#
#################################################################################### OUTPUTS ####
#
# [1] supplementary_ROSMAP/Results/GSEA/
#     GSEA results.
#
############################################################################# INITIALISATION ####

require(biomaRt)
require(dplyr)
require(data.table)
require(piano)

infile_all <- c("supplementary_ROSMAP/Results/DESeq2/allDEG_c1_vs_c4.tsv",
                "supplementary_ROSMAP/Results/DESeq2/allDEG_c2_vs_c4.tsv",
                "supplementary_ROSMAP/Results/DESeq2/allDEG_c3_vs_c4.tsv")

infile_GSEASets <- "supplementary_ROSMAP/Data/GSEA/BioMartforGSC.tsv"

outTable <- c("supplementary_ROSMAP/Results/GSEA/GSEA_c1_vs_c4.tsv",
              "supplementary_ROSMAP/Results/GSEA/GSEA_c2_vs_c4.tsv",
              "supplementary_ROSMAP/Results/GSEA/GSEA_c3_vs_c4.tsv")

outHeatmap <- c("supplementary_ROSMAP/Results/GSEA/GSEA_c1_vs_c4.png",
                "supplementary_ROSMAP/Results/GSEA/GSEA_c2_vs_c4.png",
                "supplementary_ROSMAP/Results/GSEA/GSEA_c3_vs_c4.png")

height <- c(6700, 2900, 9600)


#### Functions #############################################################################################

### Get biomaRt results
fetchBM <- function(x,                                                                                      # Character or vector to submit to biomaRt
                    file = NULL                                                                             # File to check for reading or writing
) {
    # If file is specified, see if biomaRt results already exist, and if so, load them
    if (!is.null(file)) {                                                                                   # Check if file is specified
        try(
            suppressWarnings(res <- read.delim(file = file,                                                 # Try to load the file
                                               colClasses = "character")),
            silent = TRUE
        )
    }
    # If biomaRt results don't exist locally, fetch them from biomaRt and save them to file if specified
    while (!exists("res")) {                                                                                # Check if biomaRt results were loaded successfully
        try(                                                                                                # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("ensembl_gene_id",
                                        "name_1006"),
                         filters = "ensembl_gene_id",
                         values = x,
                         mart = useMart("ensembl",
                                        dataset = "hsapiens_gene_ensembl")))
        if (exists("res") & (!is.null(file))) {                                                             # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

#### Main ##################################################################################################

doGSEA <- function(infile, outHeatmap, outTable, height) {

    data <- fread(infile)

    # Collect directions

    data <- data %>%
        filter(!is.na(pval)) %>%                                                                                # Remove na in pval
        mutate(dir = sign(l2fc))

    dir <- data$dir
    names(dir) <- data$gene

    # Collect p values

    pval <- data$pval
    names(pval) <- data$gene

    # Query biomaRt to make gene set collections

    sets <- fetchBM(x = data$gene %>% unique(),
                    file = infile_GSEASets) %>%
        filter(name_1006 != "") %>%
        loadGSC()

    # Run GSEA

    gsaRes <- runGSA(geneLevelStats = pval,
                     directions = dir,
                     gsc = sets)

    # Save tables
    fwrite(x = GSAsummaryTable(gsaRes),
           file = outTable,
           sep = "\t")

    # Plot heatmaps

    png(file = outHeatmap,
        width = 940,
        height = height)
    print({
        GSAheatmap(gsaRes,
                   ncharLabel = 64,
                   cellnote = "none",
                   colorkey = F,
                   cex = 1.5)
    })

    dev.off()
}

for (i in 1:length(infile_all)) {
    doGSEA(infile_all[[i]], outHeatmap[[i]], outTable[[i]], height[[i]])
}

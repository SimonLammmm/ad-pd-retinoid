################################################################################### METADATA ####
#
# zeb_consensusExpression.R
# Version: 1.0.0 (2020-10-20)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Determine mean expression counts per genotype.
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
# [-] zeb/Data/ConsensusExpression/bmindex.tsv
#     Optional biomaRt results dictionary. Will be created dynamically if not provided.
#
##################################################################################### INPUTS ####
#
# [1] zeb/Data/finalCts/colData.tsv
#     Metadata.
#
# [2] zeb/Data/finalCts/finalCts.tsv
#     Cleaned and combined counts ready for submission to summarizedExperiment.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Results/ConsensusExpression
#     Mean expression counts per genotype.
#
############################################################################# INITIALISATION ####

require(dplyr)
require(data.table)
require(biomaRt)

incounts <- "zeb/Data/finalCts/finalCts.tsv"
incoldata <- "zeb/Data/finalCts/colData.tsv"
outdir <- "zeb/Results/ConsensusExpression/"
bmindex <- "zeb/Data/ConsensusExpression/bmindex.tsv"

################################################################################## FUNCTIONS ####

### Subset counts across organs

subsetCts <- function(cts, colData, organSel, strainSel) {
    colData <- colData %>%
        filter(part == organSel) %>%
        filter(tert1 == strainSel)

    cts %>%
        dplyr::select(entrezgene_id, colData$sample)
}


### Calculate mean counts

meanCts <- function(subsetCts) {
    cbind(subsetCts$entrezgene_id, rowMeans(subsetCts[,2:length(subsetCts)])) %>%
        as.data.frame()
}


### biomaRt: convert ENSDART to entrez
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
        try({                                                                                                # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("entrezgene_id",
                                        "ensembl_transcript_id"),
                         filters = "ensembl_transcript_id",
                         values = x,
                         mart = useMart("ensembl",
                                        dataset = "drerio_gene_ensembl"))
            Sys.sleep(1)
            })
        if (exists("res") & (!is.null(file))) {                                                             # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}



### Input ##################################################################################################

cts <- fread(incounts)                                                                                      # Read counts file

genes_ENSDARG <- fetchBM(cts$target_id, file = bmindex) %>%                                                 # Convert ENSDART to entrez
    filter(!is.na(entrezgene_id))
cts <- cts %>%
    right_join(genes_ENSDARG, by = c("target_id" = "ensembl_transcript_id")) %>%
    dplyr::select(entrezgene_id, 2:length(cts))

colData <- fread(incoldata)                                                     # Read coldata




### Run ####################################################################################################

### Subset

cts_Full_wt <- subsetCts(cts, colData, c("Brain", "Muscle", "Skin", "Liver"), "wt")
cts_Full_het <- subsetCts(cts, colData, c("Brain", "Muscle", "Skin", "Liver"), "het")
cts_Full_dko <- subsetCts(cts, colData, c("Brain", "Muscle", "Skin", "Liver"), "dko")

cts_Brain_wt <- subsetCts(cts, colData, "Brain", "wt")
cts_Brain_het <- subsetCts(cts, colData, "Brain", "het")
cts_Brain_dko <- subsetCts(cts, colData, "Brain", "dko")

cts_Muscle_wt <- subsetCts(cts, colData, "Muscle", "wt")
cts_Muscle_het <- subsetCts(cts, colData, "Muscle", "het")
cts_Muscle_dko <- subsetCts(cts, colData, "Muscle", "dko")

cts_Skin_wt <- subsetCts(cts, colData, "Skin", "wt")
cts_Skin_het <- subsetCts(cts, colData, "Skin", "het")
cts_Skin_dko <- subsetCts(cts, colData, "Skin", "dko")

cts_Liver_wt <- subsetCts(cts, colData, "Liver", "wt")
cts_Liver_het <- subsetCts(cts, colData, "Liver", "het")
cts_Liver_dko <- subsetCts(cts, colData, "Liver", "dko")


### Calculate mean counts

meanCts_Full_wt <- meanCts(cts_Full_wt)
meanCts_Full_het <- meanCts(cts_Full_het)
meanCts_Full_dko <- meanCts(cts_Full_dko)

meanCts_Brain_wt <- meanCts(cts_Brain_wt)
meanCts_Brain_het <- meanCts(cts_Brain_het)
meanCts_Brain_dko <- meanCts(cts_Brain_dko)

meanCts_Muscle_wt <- meanCts(cts_Muscle_wt)
meanCts_Muscle_het <- meanCts(cts_Muscle_het)
meanCts_Muscle_dko <- meanCts(cts_Muscle_dko)

meanCts_Skin_wt <- meanCts(cts_Skin_wt)
meanCts_Skin_het <- meanCts(cts_Skin_het)
meanCts_Skin_dko <- meanCts(cts_Skin_dko)

meanCts_Liver_wt <- meanCts(cts_Liver_wt)
meanCts_Liver_het <- meanCts(cts_Liver_het)
meanCts_Liver_dko <- meanCts(cts_Liver_dko)


### Write to file

fwrite(meanCts_Full_wt, paste0(outdir, "meanCts_full_wt.tsv"), sep = "\t")
fwrite(meanCts_Full_het, paste0(outdir, "meanCts_full_het.tsv"), sep = "\t")
fwrite(meanCts_Full_dko, paste0(outdir, "meanCts_full_dko.tsv"), sep = "\t")

fwrite(meanCts_Brain_wt, paste0(outdir, "meanCts_brain_wt.tsv"), sep = "\t")
fwrite(meanCts_Brain_het, paste0(outdir, "meanCts_brain_het.tsv"), sep = "\t")
fwrite(meanCts_Brain_dko, paste0(outdir, "meanCts_brain_dko.tsv"), sep = "\t")

fwrite(meanCts_Muscle_wt, paste0(outdir, "meanCts_muscle_wt.tsv"), sep = "\t")
fwrite(meanCts_Muscle_het, paste0(outdir, "meanCts_muscle_het.tsv"), sep = "\t")
fwrite(meanCts_Muscle_dko, paste0(outdir, "meanCts_muscle_dko.tsv"), sep = "\t")

fwrite(meanCts_Skin_wt, paste0(outdir, "meanCts_skin_wt.tsv"), sep = "\t")
fwrite(meanCts_Skin_het, paste0(outdir, "meanCts_skin_het.tsv"), sep = "\t")
fwrite(meanCts_Skin_dko, paste0(outdir, "meanCts_skin_dko.tsv"), sep = "\t")

fwrite(meanCts_Liver_wt, paste0(outdir, "meanCts_liver_wt.tsv"), sep = "\t")
fwrite(meanCts_Liver_het, paste0(outdir, "meanCts_liver_het.tsv"), sep = "\t")
fwrite(meanCts_Liver_dko, paste0(outdir, "meanCts_liver_dko.tsv"), sep = "\t")

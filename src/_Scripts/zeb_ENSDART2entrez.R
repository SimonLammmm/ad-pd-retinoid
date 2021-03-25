################################################################################### METADATA ####
#
# zeb_ENSDART2entrez.R
# Version: 1.0.0 (2020-10-22)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Convert ensembl gene symbols to entrez gene symbols in DESeq2 results.
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
# [1] zeb/Results/DESeq2/sep/
#     DESeq2 results.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Results/DESeq2/sep/entrez/
#     DESeq2 results.
#
############################################################################# INITIALISATION ####

#### Attach packages ####
require(biomaRt)
require(dplyr)
require(data.table)


#### User variables ####
indir <- "zeb/Results/DESeq2/all/sep/"
outdir <- "zeb/Results/DESeq2/all/sep/entrez/"
bmindex <- "zeb/Data/ConsensusExpression/bmindex.tsv"


################################################################################## FUNCTIONS ####

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

bmdata <- fetchBM(data$gene, bmindex)

ENSDART2entrez <- function(infile, outfile) {
    data <- fread(infile)
    data_entrez <- left_join(data, bmdata, by = c("gene" = "ensembl_transcript_id"))
    data_entrez <- dplyr::select(data_entrez, gene = entrezgene_id, 2:length(data_entrez))
    fwrite(data_entrez, outfile, sep = "\t")
}

####################################################################################### MAIN ####

infile <- list.files(indir, pattern = "DESeq2")
outfile <- paste0(outdir, infile)

for (i in 1:length(infile)) {
    try({
        ENSDART2entrez(paste0(indir, infile[i]), outfile[i])
    })
}

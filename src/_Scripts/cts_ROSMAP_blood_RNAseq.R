################################################################################### METADATA ####
#
# cts_ROSMAP_blood_RNAseq.R
# Version: 1.0.1 (2020-08-13)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Conversion  of  kallisto  counts  file  to   tab-separated   values  counts  files  for  ROSMAP
# blood RNA-seq.
#
# Citation:
# Bennett DA  et al, 2012,  "Overview  and  findings  from the  religious  orders  study."  Curr.
# Alzheimer Res., 9(6): 628-645. DOI: 10.2174/156720512801322573.
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
##################################################################################### INPUTS ####
#
# [1] cts_ROSMAP_blood_RNAseq/Data/Sample_XXX_END1.fastq.gz-trimmed.fastq.gz-kallisto_out/
#     Directories containing kallisto quant results
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_blood_RNAseq/Results/
#     Expression counts in tsv format.
#
############################################################################# INITIALISATION ####

require(foreach)
require(dplyr)
require(data.table)
require(biomaRt)

indir <- "cts_ROSMAP_blood_RNAseq/Data/"
outdir <- "cts_ROSMAP_blood_RNAseq/Results/"

################################################################################## FUNCTIONS ####

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
            res <- getBM(attributes = c("ensembl_gene_id",
                                        "ensembl_transcript_id"),
                         filters = "ensembl_transcript_id",
                         values = x,
                         mart = useMart("ensembl",
                                        dataset = "hsapiens_gene_ensembl")))
        if (exists("res") & (!is.null(file))) {                                      # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

####################################################################################### MAIN ####

infiles <- paste0(list.files(indir, full.names = T, pattern = "kallisto_out"),                                                        # load estimated abundance files
                  "/abundance.tsv")

infiles_samplenames <- list.files(indir) %>%                                                                # record friendly sample names
    sub(pattern = "_E.+", replacement = "")

data <- foreach (file = infiles,                                                                            # record target id (version nunber removed) and tpm
                 .final = function(file) setNames(file, infiles_samplenames)) %do% {
                     fread(file) %>%
                         dplyr::transmute(target_id = sub(x = target_id, pattern = "\\.\\d+", replacement = ""),
                                          tpm)
                 }

genes_ENST <- foreach (s = data) %do% {                                                                     # get the list of genes (ENST)
    s %>%
        dplyr::select(target_id)
} %>%
    unlist() %>%
    unique()

genes_ENSG <- fetchBM(genes_ENST)                                                                           # determine the list of genes (ENSG)

data2 <- foreach (s = names(data),                                                                          # join ENSTs, remove target_id column
                  .final = function(s) setNames(s, names(data))) %do% {
                      d <- data[[s]] %>%
                          inner_join(y = genes_ENSG,
                                     by = c("target_id" = "ensembl_transcript_id")) %>%
                          dplyr::select(ensembl_gene_id,
                                    tpm)
                      names(d)[2] <- s
                      d
                  }

for (s in names(data2)) {                                                                                   # write to file
    fwrite(x = data2[[s]],
           file = paste0(outdir,
                         s,
                         ".tsv"),
           sep = "\t")
}

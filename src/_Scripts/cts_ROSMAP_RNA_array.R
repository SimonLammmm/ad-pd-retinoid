################################################################################### METADATA ####
#
# cts_ROSMAP_RNA_array.R
# Version: 1.0.1 (2020-07-27)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Parsing of ROSMAP RNA array to tab-separated values format.
#
# Citation:
# Bennett DA  et al, 2012,  "Overview  and  findings  from the  religious  orders  study."  Curr.
# Alzheimer Res., 9(6): 628-645. DOI: 10.2174/156720512801322573.
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
# [1] cts_ROSMAP_RNA_array/Data/ROSMAP_arrayExpression_normalized.tsv
#     https://www.synapse.org/#!Synapse:syn4009614
#     Approval required.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_RNA_array/Results/
#     Expression counts in tsv format.
#
############################################################################# INITIALISATION ####
require(foreach)
require(dplyr)
require(data.table)
require(biomaRt)

infile <- "cts_ROSMAP_RNA_array/Data/ROSMAP_arrayExpression_normalized.tsv"
outdir <- "cts_ROSMAP_RNA_array/Results/"

####################################################################################### MAIN ####

dir.create(outdir, showWarnings = F)

#####

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
                                        "illumina_humanht_12_v4"),
                         filters = "illumina_humanht_12_v4",
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

#####


data <- fread(infile)

data <- data %>%
    dplyr::select(-ProbeID, -Symbol)

genes <- fetchBM(data$TargetID)

data <- data %>%
    inner_join(y = genes,
              by = c("TargetID" = "illumina_humanht_12_v4")) %>%
    dplyr::select(-TargetID)

data2 <- foreach (s = names(data),
                  .final = function(s) setNames(s, names(data))) %do%
                  {
                      if (s == "ensembl_gene_id") {
                          NULL
                      }
                      else {
                          data %>%
                              dplyr::select(ensembl_gene_id,
                                     s)
                      }
                  }

data2$ensembl_gene_id <- NULL

for (s in names(data2)) {
    fwrite(x = data2[[s]],
           file = paste0(outdir,
                         s,
                         ".tsv"),
           sep = "\t")
}

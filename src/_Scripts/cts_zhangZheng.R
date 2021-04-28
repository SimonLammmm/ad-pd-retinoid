################################################################################### METADATA ####
#
# cts_zhangZheng.R
# Version: 1.0.0 (2020-03-03)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Conversion of GEO GSE20295 Affymetrix microarray data to ensembl counts.
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
########################################################################## FILE REQUIREMENTS ####
#
#
# [1] cts_zhangZheng/Data/GSE20295_RAW/
#     https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE20295&format=file
#
# [-] biomaRt results: GEM genes (default = cts_zhangZheng/Data/biomaRt_affyMappings.tsv) [optional file]
#     Optional file. Will be dynamically created if absent.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_zhangZheng/Results/
#     Expression counts in tsv format.
#
############################################################################# INITIALISATION ####

require(affy)
require(biomaRt)
require(data.table)
require(doParallel)
require(dplyr)
require(foreach)
require(hgu133acdf)

registerDoParallel(detectCores() - 2)

userFile_affyMappings <- "cts_zhangZheng/Data/biomaRt_affyMappings.tsv"
userDir_affyCEL <- "cts_zhangZheng/Data/GSE20295_RAW/"

outdir <- "cts_zhangZheng/Results/"

################################################################################## FUNCTIONS ####

### Start biomaRt
initBM <- function() {
    hosts <- c(                                                                      # Comment out lines if you don't want to use a particular mirror
        "www.ensembl.org"
        ,
        "useast.ensembl.org"
        ,
        "uswest.ensembl.org"
        ,
        "asia.ensembl.org"
    )
    useMart("ensembl",
            dataset = "hsapiens_gene_ensembl",
            host = hosts[sample(1:length(hosts),                                     # Choose a mirror at random (out of the ones that weren't commented out)
                                1)])
}

### Get biomaRt results (from ensembl gene id)
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
        try(
            mart <- initBM()
        )
        try(                                                                         # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("ensembl_gene_id",
                                        "hgnc_symbol",
                                        "affy_hg_u133a",
                                        "description"),
                         filters = "affy_hg_u133a",
                         values = x,
                         mart = mart))
        if (exists("res") & (!is.null(file))) {                                      # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

####################################################################################### MAIN ####


g <- fetchBM(x = names(hgu133acdf),
             file = userFile_affyMappings)

expr <- ReadAffy(filenames = list.files(path = userDir_affyCEL,
                                        pattern = "*.CEL$",
                                        full.names = T)) %>%
    rma()

expr <- tibble(probe = rownames(expr@assayData$exprs)) %>%
    cbind(expr@assayData$exprs) %>%
    as_tibble() %>%
    inner_join(y = g,
               by = c("probe" = "affy_hg_u133a"))

expr <- foreach(a = grep(pattern = "\\.CEL$",
                         x = names(expr)),
                .final = function(a) setNames(a, regmatches(x = names(expr),
                                                            m = regexpr(pattern = "^GSM\\d{6}",
                                                                        text = names(expr))))) %dopar%
                                                                        {
                                                                            tibble(gene = expr[["ensembl_gene_id"]],
                                                                                   expr = expr[[a]])
                                                                        }


dir.create(outdir,
           showWarnings = F)

foreach(a = names(expr)) %dopar% {
    fwrite(x = expr[[a]],
           file = paste(outdir,
                        a,
                        "_expr.tsv",
                        sep = ""),
           sep = "\t")
}

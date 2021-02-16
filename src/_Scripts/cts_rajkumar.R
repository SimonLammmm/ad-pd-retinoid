################################################################################### METADATA ####
#
# cts_rajkumar.R
# Version: 1.0.0 (2020-05-06)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Conversion of gene symbol counts in Microsoft Excel counts to ensembl gene counts in  tab-
# separated values format for the processing of data reported by Rajkumar et al.
#
# Citation:
# Rajkumar AP  et al, 2020, "Postmortem  cortical  transcriptomics  of Lewy body  dementia reveal
# mitochondrial dysfunction  and lack of neuroinflammation."  Am. J. Geriatr.  Psychiatry, 28(1):
# 75-86. DOI: 10.1016/j.jagp.2019.06.007.
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
#
# [1] cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc2.xlsx
#     https://ars.els-cdn.com/content/image/1-s2.0-S1064748119304142-mmc2.xlsx
#
# [2] cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc3.xlsx
#     https://ars.els-cdn.com/content/image/1-s2.0-S1064748119304142-mmc3.xlsx
#
#
#################################################################################### OUTPUTS ####
#
# [1] cts_rajkumar/Results/
#     Expression counts in tsv format.
#
############################################################################# INITIALISATION ####

require(readxl)
require(dplyr)
require(foreach)
require(data.table)
require(biomaRt)
require(tidyr)
require(tidyselect)

userPath_acc <- "cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc2.xlsx"
userPath_dlpfc <- "cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc3.xlsx"

outdir <- "cts_rajkumar/Results/"

################################################################################## FUNCTIONS ####


getBioMart <- function(genes, filters) {                                  # Convert symbol to geneid with biomaRt
    getBM(attributes = c("ensembl_gene_id",
                         filters),
          filters = filters,
          values = genes,
          mart = useMart("ensembl",
                         dataset = "hsapiens_gene_ensembl",
                         host = "useast.ensembl.org"))              # Querying human genes

}

####################################################################################### MAIN ####

cts_acc <- read_xlsx(path = userPath_acc)               # Read in anterior cingulate cortex counts
cts_dlpfc <- read_xlsx(path = userPath_dlpfc)           # Read in dorsolateral prefrontal cortex counts

names(cts_acc)[-1] <- names(cts_acc)[-1] %>%            # Append _acc to sample names
    sub(pattern = "$",
        replacement = "_acc")

names(cts_dlpfc)[-1] <- names(cts_dlpfc)[-1] %>%        # Append _dlpfc to sample names
    sub(pattern = "$",
        replacement = "_dlpfc")

cts <- full_join(x = cts_acc,                           # Merge lists
                 y = cts_dlpfc,
                 by = "Gene")

filters <- c("clone_based_ensembl_gene",              # Convert gene names in list to ENSG: define name conventions used in the gene list
             "hgnc_symbol",
             "entrezgene_accession",
             "external_synonym")


geneids <- foreach (f = filters,                        # Query biomaRt
                    .final = function(f) setNames(f, filters)) %do% {
                        getBioMart(genes = cts$Gene,
                                   filters = f)
                    }

genes <- tibble(ensembl_gene_id = "")[0,]               # Collapse list to ENSG
for (f in filters) {
    genes <- genes %>%
        full_join(y = geneids[[f]],
                  by = "ensembl_gene_id")
}

genes <- genes %>%
    gather(cols = filters,
           key = "Source",
           value = "Gene")

cts <- cts %>%                                          # Merge ENSG to counts
    left_join(y = genes,
              by = "Gene") %>%
    dplyr::select(ensembl_gene_id, 2:41) %>%
    unique() %>%                                        # Remove duplicates
    filter(ensembl_gene_id != "")                                # Remove empties


dir.create(outdir,
           showWarnings = F)

foreach (n = names(cts)) %do% {
    if (n != "ensembl_gene_id") {
        fwrite(x = cts[c("ensembl_gene_id", n)],
               file = paste(outdir,
                            n,
                            "_expr.tsv",
                            sep = ""),
               sep = "\t")
    }
}

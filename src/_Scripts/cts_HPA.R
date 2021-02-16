################################################################################### METADATA ####
#
# cts_HPA.R
# Version: 1.0.1 (2020-05-18)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Conversion and compartmentation of HPA counts to  brain regions of three levels of granularity:
# broad, precise, and single-cell.
#
# Citation:
# Uhlén M et al, 2015, "Tissue-based map of the human proteome." Science 347(6220): 1260419. DOI:
# 10.1126/science.1260419.
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
# [1] cts_HPA/Data/normal_tissue.tsv
#     https://www.proteinatlas.org/download/normal_tissue.tsv.zip
#
#################################################################################### OUTPUTS ####
#
# [1] cts_HPA/Results/broad
# [2] cts_HPA/Results/precise
# [3] cts_HPA/Results/sc
#     Expression counts in tsv format.
#
############################################################################# INITIALISATION ####

require(data.table)
require(dplyr)
require(foreach)

userFile_In <- "cts_HPA/Data/normal_tissue.tsv"

outdir_broad <- "cts_HPA/Results/broad/"
outdir_precise <- "cts_HPA/Results/precise/"
outdir_sc <- "cts_HPA/Results/sc/"

####################################################################################### MAIN ####

#### Counts matrix

tissues <- c("caudate",              # These brain tissues are covered in HPA normal_tissue.tsv
             "cerebellum",
             "cerebral cortex",
             "hippocampus",
             "hypothalamus",
             "pituitary gland",
             "dorsal raphe",
             "choroid plexus",
             "substantia nigra")

cts <- fread(userFile_In) %>%                              # Read in the HPA data
    filter(Tissue %in% tissues) %>%                        # Return only brain tissue data
    mutate(expr = case_when(Level == "Not detected" ~ 0,   # Re-level staining strengths to numerical data
                            Level == "Low" ~ 1,
                            Level == "Medium" ~ 2,
                            Level == "High" ~ 3),
           broadTissue = case_when(Tissue == "caudate" | Tissue == "substantia nigra" ~ "basal ganglia",
                                   Tissue == "cerebellum" ~ "cerebellum",
                                   Tissue == "cerebral cortex" ~ "cerebral cortex",
                                   Tissue == "hippocampus" ~ "hippocampus",
                                   Tissue == "hypothalamus" | Tissue == "pituitary gland" ~ "hypothalamus",
                                   Tissue == "dorsal raphe" ~ "pons and medulla",
                                   Tissue == "choroid plexus" ~ "choroid plexus"))

getList <- function(x) {
    x %>% unique() %>% sort()
}

cts_tissueB <- getList(cts$broadTissue)                # Get the list of all broad regions in the data
cts_tissue <- getList(cts$Tissue)                      # Get the list of all precise regions in the data
cts_celltype <- getList(cts$`Cell type`)               # Get the list of all single cells in the data

#### Single-cell models

sc <- foreach (tis = cts_tissue) %:%              # Compile single cell models for each tissue
    foreach (cel = cts_celltype,
             .final = function(cel) setNames(cel, paste0(tis, "_", cts_celltype))) %do% {
        cts %>%
            filter(Tissue == tis) %>%
            filter(`Cell type` == cel)
    } %>%
    unlist(recursive = F)

sc[lapply(X = sc,
          FUN = function(x) nrow(x) == 0) %>%   # Remove empties
       as.logical()] <- NULL

#### Precise models

prec <- foreach (tis = cts_tissue,
                 .final = function(tis) setNames(tis, cts_tissue)) %do% {
                     cts %>%
                         filter(Tissue == tis)
                 }


#### Broad models

broad <- foreach (tisb = cts_tissueB,
                  .final = function(tisb) setNames(tisb, cts_tissueB)) %do% {
                      cts %>%
                          filter(broadTissue == tisb)
                  }


#### Output

dir.create(outdir_broad)
dir.create(outdir_precise)
dir.create(outdir_sc)

foreach (m = names(sc)) %do% {
    fwrite(x = sc[[m]],
           file = paste0(outdir_sc,
                         m,
                         ".tsv"),
           sep = "\t")
}

foreach (m = names(prec)) %do% {
    fwrite(x = prec[[m]],
           file = paste0(outdir_precise,
                         m,
                         ".tsv"),
           sep = "\t")
}

foreach (m = names(broad)) %do% {
    fwrite(x = broad[[m]],
           file = paste0(outdir_broad,
                         m,
                         ".tsv"),
           sep = "\t")
}

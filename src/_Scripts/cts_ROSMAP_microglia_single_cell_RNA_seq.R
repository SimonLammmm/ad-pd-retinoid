################################################################################### METADATA ####
#
# cts_ROSMAP_microglia_single_cell_RNA_seq.R
# Version: 1.0.0 (2020-08-03)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Conversion  of cellranger  counts  file  to   tab-separated  values  counts  files  for  ROSMAP
# microglia single-cell RNA-seq.
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
# [1] cts_ROSMAP_microglia_single_cell_RNA_seq/Results/CellRanger/
#     Directories containing Cell Ranger counts
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_microglia_single_cell_RNA_seq/Results/Expression/
#     Expression counts in tsv format
#
############################################################################# INITIALISATION ####
require(R.utils)
require(data.table)
require(Matrix)
require(dplyr)
require(foreach)

indir <- "cts_ROSMAP_microglia_single_cell_RNA_seq/Results/CellRanger/"
outdir <- "cts_ROSMAP_microglia_single_cell_RNA_seq/Results/Expression/"

####################################################################################### MAIN ####

dir.create(outdir, showWarnings = F)

indirs <- list.dirs(indir, recursive = F)

indirs_names <- indirs %>%
    sub(pattern = ".+(Microglia_MO_.+-counts).*",
        replacement = "\\1")

inpath_barcode <- paste0(indirs, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
inpath_features <- paste0(indirs, "/outs/filtered_feature_bc_matrix/features.tsv.gz")
inpath_matrix <- paste0(indirs, "/outs/filtered_feature_bc_matrix/matrix.mtx.gz")



data <- foreach(f = 1:length(indirs),                                                                       # import to sparse matrix
                .final = function(f) setNames(f, indirs_names)) %do% {
                    m <- readMM(file = inpath_matrix[f])
                    r <- fread(inpath_features[f], header = F)
                    c <- fread(inpath_barcode[f], header = F)
                    rownames(m) <- r[[1]]
                    colnames(m) <- c[[1]]
                    m
                }

consensusData <- foreach(f = names(data),                                                                   # calculate means for each gene
                         .final = function(f) setNames(f, names(data))) %do% {
                             d <- data.frame(row.names(data[[f]]), rowMeans(data[[f]]))
                             names(d) <- c("gene", f)
                             d
                         }

foreach(f = names(consensusData)) %do% {
    fwrite(x = consensusData[[f]],
           file = paste0(outdir, f, ".tsv"),
           sep = "\t")
}

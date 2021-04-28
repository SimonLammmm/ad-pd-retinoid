################################################################################### METADATA ####
#
# zeb_DESeq2_sep.R
# Version: 1.0.0 (2020-11-06)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Separate DESeq2 results files into separate files.
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
# [1] zeb/Results/DESeq2
#     DESeq2 results.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Results/DESeq2
#     DESeq2 results.
#
############################################################################# INITIALISATION ####
require(dplyr)
require(data.table)
require(foreach)

inDESeq1 <- "zeb/Results/DESeq2/all/zeballGenes_ByGenotype.tsv"
inDESeq2 <- "zeb/Results/DESeq2/all/zeballGenes_ByPart.tsv"

outdir <-  "zeb/Results/DESeq2/all/sep/"

data1 <- fread(inDESeq1)
data2 <- fread(inDESeq2)

################################################################################## FUNCTIONS ####

subsetter <- function(.data) {

    foreach (i = unique(.data$subset),
             .final = function(i) setNames(i, unique(.data$subset))) %do% {
                 foreach (j = unique(.data$differential),
                          .final = function(j) setNames(j, unique(.data$differential))) %do% {
                              foreach (k = unique(.data$reference),
                                       .final = function(k) setNames(k, unique(.data$reference))) %do% {
                                           if (j != k) {
                                               .data %>%
                                                   filter(subset == i) %>%
                                                   filter(differential == j) %>%
                                                   filter(reference == k)
                                           }
                                       }
                          }

             }
}

saver <- function(.data, outdir) {
    for (i in 1:length(.data)) {
        for (j in 1:length(.data[[i]])) {
            for (k in 1:length(.data[[i]][[j]])) {
                if (names(.data[[i]])[j] != names(.data[[i]][[j]])[k]) {
                    fwrite(.data[[i]][[j]][[k]],
                           paste0(outdir,
                                  "DESeq2_",
                                  names(.data)[i],
                                  "_",
                                  names(.data[[i]])[j],
                                  "_vs_",
                                  names(.data[[i]][[j]])[k],
                                  ".tsv"),
                           sep = "\t")
                }
            }
        }
    }
}

####################################################################################### MAIN ####

data1 <- subsetter(data1)
data2 <- subsetter(data2)

saver(data1, outdir)
saver(data2, outdir)

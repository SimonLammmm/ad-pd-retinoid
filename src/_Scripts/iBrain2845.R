################################################################################### METADATA ####
#
# iBrain2845.R
# Version: 2.2.10 (2020-01-22)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Generate a generalised  genome-scale metabolic model (GEM)  based on a context-specific GEM and
# a general database.  Use the reactions  from the context-specific  GEM and generalise the model
# using the genes from HMR3 that catalyse those reactions.
#
# Citation:
# Mardinoglu A et al, 2013, "Integration of clinical  data with a genome-scale metabolic model of
# the human adipocyte." Mol. Syst. Biol. 9: 649.
# Mardinoglu A  et al, 2014,  "Genome-scale  metabolic  modelling of  hepatocytes  reveals serine
# deficiency in patients with non-alcoholic fatty liver disease." Nat. Commun. 5: 3083.
#
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
# Modify this script with  the filepaths of the below required and  optional files in the section
# entitled "User variables".  Optional user variables  (i.e. those with 'biomaRt'  in their name)
# may be omitted by setting these variables to NULL.
#
# [1] iBrain2845/Data/iAdipocytes1850_v2.xlsx
#     Manually curated from Mardinoglu et al, 2013.
#
# [2] iBrain2845/Data/HMRdatabase3_00.xlsx
#     Manually curated from Mardinoglu et al, 2014.
#
#################################################################################### OUTPUTS ####
#
# [1] iBrain2845/Results/iBrain2845.sbml.xlsx
#     Generalised GEM.
#
################################################################################### OVERVIEW ####
#
# 0. INITIALISE ENVIRONMENT
#    a. Attach packages
#    b. Define functions
#    c. Set user variables
# 1. RXNS sheet
# 2. METS sheet
# 3. COMPS sheet
# 4. MODELS sheet
# 5. GENES sheet
# 6. ARTIFICIAL sheet
# 7. Write to Excel
#
############################################################################# INITIALISATION ####

### Attach required packages
suppressPackageStartupMessages({
    require(readxl)                                                             # For reading Excel documents
    require(writexl)                                                            # For writing to Excel format
})


############################################################################# USER VARIABLES ####

userFile_generalGEM <- "iBrain2845/Data/HMRdatabase3_00.xlsx"                                   # Path to the general GEM file [1]
userFile_contextGEM <- "iBrain2845/Data/iAdipocytes1850_v2-withRxnNames.xlsx"                 # Path to the context-specific GEM file [2]
userFile_outputGEM <- "iBrain2845/Results/iBrain2845.sbml.xlsx"                           # Directory where the output GEM should be saved


####################################################################################### MAIN ####

context <- list()
general <- list()
output <- list()

## 1. Sheet RXNS
general$RXNS <- read_excel(path = userFile_generalGEM,                          # Read the general GEM RXNS sheet
                           sheet = "RXNS",
                           col_types = "text")
context$RXNS <- read_excel(path = userFile_contextGEM,                          # Read the context-specific GEM RXNS sheet
                           sheet = "RXNS",
                           col_types = "text")
output$RXNS <- context$RXNS                                                     # Accept all reactions in the context-specific GEM
output$RXNS$`#`[which(!(context$RXNS$RXNID %in% general$RXNS$RXNID))] <- "#"    # Comment out reactions that are only in the context-specific GEM
output$RXNS$`GENE ASSOCIATION`[which(context$RXNS$RXNID %in% general$RXNS$RXNID)] <-   # Replace the gene associations of all context-specific GEM reactions with those from the general GEM
    general$RXNS$`GENE ASSOCIATION`[which(general$RXNS$RXNID %in% context$RXNS$RXNID)]

## 2. Sheet METS
output$METS <- read_excel(path = userFile_contextGEM,                           # Read the context-specific GEM METS sheet and designate it to output GEM
                          sheet = "METS",
                          col_types = "text")
mets_unlist <- unique(sub(x = grep(x = unique(unlist(strsplit(x = output$RXNS$EQUATION[which(is.na(output$RXNS$`#`))], # Gather the list of metabolites actually used in output GEM reactions without "#"
                                                              split = "\\s(\\+|=>|<=>)\\s"))),  # Split reaction strings by whitespace, then literal + or => or <=>, then whitespace to get metabolites with stoichiometries
                                   pattern = "\\[.\\]$",
                                   value = TRUE),
                          pattern = "^(\\d+|\\d+\\.\\d+)\\s",                   # Remove integers and decimals at the start of the string to remove stoichiometries from metabolite names
                          replacement = ""))

output$METS$`#`[!(output$METS$METID %in% mets_unlist)] <- "#"                   # Fill the comment column of the unused metabolites in output GEM with "#"

## 3. Sheet COMPS
output$COMPS <- read_excel(path = userFile_contextGEM,                          # Read the context-specific GEM COMPS sheet and designate it to output GEM
                           sheet = "COMPS",
                           col_types = "text")

## 4. Sheet MODEL
modelName_output <- paste(userFile_outputGEM,                                   # Set the output GEM full model name
                          as.character(nrow(output$RXNS) - nrow(output$DELRXNS)),
                          sep = "")
output$MODEL <- data.frame(`#` = "",                                            # Designate model metadata to output GEM
                           MODELID = modelName_output,
                           MODELNAME = "Genome-scale metabolic model for brain",
                           `DEFAULT LOWER` = -1000,
                           `DEFAULT UPPER` = 1000,
                           `CONTACT GIVEN NAME` = "Simon",
                           `CONTACT FAMILY NAME` = "Lam",
                           `CONTACT EMAIL` = "simon.1.lam@kcl.ac.uk",
                           `ORGANIZATION` = "King's College London",
                           `TAXONOMY` = "taxonomy:9606",
                           `NOTES` = '<p>This SBML representation of the Homo sapiens brain metabolic network is made available under the Creative Commons Attribution-Share Alike 3.0 Unported Licence (see <a href="http://www.creativecommons.org">www.creativecommons.org</a>).</p>',
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

## 5. Sheet GENES
general$GENES <- read_excel(path = userFile_generalGEM,                         # Read the general GEM GENES sheet
                            sheet = "GENES",
                            col_types = "text")
context$GENES <- read_excel(path = userFile_contextGEM,                         # Read the context-specific GEM GENES sheet
                            sheet = "GENES",
                            col_types = "text")
genes_unlist <- unique(unlist(strsplit(x = output$RXNS$`GENE ASSOCIATION`[which(is.na(output$RXNS$`#`))],
                                       split = ";")))                           # Store the list of genes in reactions of output GEM that are not "#"
output$GENES <- general$GENES[general$GENES$`GENE NAME` %in% genes_unlist, ]    # Accept the genes in the list from general GEM and designate to output GEM
delgenes_unlist <- !(context$GENES$`GENE NAME` %in% genes_unlist)               # Store the list of genes used in context-specific GEM but not output GEM
output$DELGENES <- context$GENES[delgenes_unlist, ]                             # Designate a list of unused genes to output GEM
try({
    output$DELGENES$`#` <- "#"                                                  # If there are any, fill the comment column of unused genes in output GEM with "#"
    output$GENES <- rbind(output$GENES,                                         # Append unused genes list to used genes list in output GEM
                          cbind(output$DELGENES,
                                data.frame(`NOTES ABOUT THE UPDATES` = "",
                                           check.names = FALSE)))
},
silent = TRUE)
output$DELGENES <- NULL                                                         # Remove the now-redundant list of deleted genes



## 6. Sheet ARTIFICIAL
output$ARTIFICIAL <- read_excel(path = userFile_contextGEM,                     # Read the context-specific GEM ARTIFICIAL sheet and designate it output GEM
                                sheet = "ARTIFICIAL",
                                col_types = "text")

## 7. Save to Excel
output <- lapply(output, as.data.frame)
write_xlsx(x = output,
           path = paste(userFile_outputGEM,
                        sep = ""))

################################################################################### METADATA ####
#
# cts_controlInteractomes-broad.R
# Version: 1.5.0 (2020-03-12)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Integration  of protein-protein  interaction  networks (PPINs)  and transcriptional  regulatory
# networks  (TRNs) as  integrated  networks  (INs) for  integration with  genome-scale  metabolic
# models (GEMs).
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
# [1] cts_controlInteractomes/Data/HuRI.tsv
#     http://www.interactome-atlas.org/data/HuRI.tsv
#
# [2] cts_controlInteractomes/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
#     https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
#
# [3] cts_controlInteractomes/Data/RegionMapping.txt
#     Tab-separated values file containing experiment  tissue label and standardised brain region
#     mappings for GTEx and FANTOM5 experiments. Manually curated.
#
# [4] cts_controlInteractomes/Data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt
#     https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz
#
#     % head cts_controlInteractomes/Data/RegionMapping.txt
#     label	preciseRegion	broadRegion
#     Brain - Amygdala	Amygdala	Amygdala
#     Brain - Anterior cingulate cortex (BA24)	Anterior cingulate cortex	Cerebral cortex
#     Brain - Caudate (basal ganglia)	Caudate	Basal ganglia
#     Brain - Cerebellar Hemisphere	Cerebellum	Cerebellum
#     Brain - Cerebellum	Cerebellum	Cerebellum
#     Brain - Cortex	Cerebral cortex	Cerebral cortex
#     Brain - Frontal Cortex (BA9)	Frontal cortex	Cerebral cortex
#     Brain - Hippocampus	Hippocampal formation	Hippocampal formation
#     Brain - Hypothalamus	Hypothalamus	Hypothalamus
#
# [4] cts_controlInteractomes/Data/FANTOM5_individual_networks/brain/
#     http://www2.unil.ch/cbg/regulatorycircuits/FANTOM5_individual_networks.tar
#
# [-] biomaRt results: GTEx genes (default = cts_controlInteractomes/Data/biomaRt_GTExTPMs.tsv) [optional file]
# [-] biomaRt restuls: HuRI genes (default = cts_controlInteractomes/Data/biomaRt_HuRI.tsv) [optional file]
#     Optional file. Will be dynamically created if absent.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_controlInteractomes/Results/broadRegion/
#     Expression counts in tsv format.
#
############################################################################ INITIALISATION ####

### Initialisation of packages

suppressPackageStartupMessages({
    require(biomaRt)               # Required for bioMart querying (initBM, fetchBM)
    require(doParallel)            # Required for parallelisation environment
    require(foreach)               # Required for parallelisation environment
    require(dplyr)                 # For data cleaning
    require(limma)                 # For normalisation of TPMs
    require(readxl)                # For reading the GEM from XLSX
})


### Set up parallelisation environment

registerDoParallel(detectCores() - 2)



################################################################################## FUNCTIONS ####

### Initialise a progress bar
pbinit <- function(x                                                                                        # Maximum progress bar value
) {
    txtProgressBar(min = 0,
                   max = x,
                   initial = 0,
                   style = 3)
}

### Start biomaRt
initBM <- function() {
    hosts <- c(                                                                                             # Comment out lines if you don't want to use a particular mirror
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
            host = hosts[sample(1:length(hosts),                                                            # Choose a mirror at random (out of the ones that weren't commented out)
                                1)])
}

### Get biomaRt results (from ensembl gene id)
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
        try(
            mart <- initBM()
        )
        try(                                                                                                # If no results were loaded successfuly, keep querying biomaRt until successful
            res <- getBM(attributes = c("ensembl_gene_id",
                                        "hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = x,
                         mart = mart))
        if (exists("res") & (!is.null(file))) {                                                             # If results were obtained from bioMart (not locally), save to file
            write.table(x = res,
                        file = file,
                        sep = "\t")
        }
    }
    res
}

############################################################################# USER VARIABLES ####

### Set user variables
userFile_HuRI             <- "cts_controlInteractomes/Data/HuRI.tsv"                                                                # Where is the HuRI master PPIN?                                              (For PPINs)
userFile_biomaRt_HuRI     <- "cts_controlInteractomes/Data/biomaRt_HuRI.tsv"                                                        # Where do you want to read/write the biomaRt mappings for HuRI from/to?      (Can be NULL)
userFile_GTExTPMs         <- "cts_controlInteractomes/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"           # Where are the GTEx TPMs?                                                    (For expression testing and PPINs)
userFile_FANTOM5TPMs      <- "cts_controlInteractomes/Data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt"                       # Where are the FANTOM5 TPMs?
userFile_RegionMapping    <- "cts_controlInteractomes/Data/RegionMapping.txt"                                                       # Where is the manually curated region mapping dictionary?                    (For setting PPIN and TRN regions)
userFile_biomaRt_GTExTPMs <- "cts_controlInteractomes/Data/biomaRt_GTExTPMs.tsv"                                                    # Where do you want to read/write the biomaRt mappings for GTEx TPMs from/to? (Can be NULL)
userDir_FANTOM5           <- "cts_controlInteractomes/Data/FANTOM5_individual_networks/brain/"                                      # Where should I look for the region label-specific FANTOM5 GGIs?             (Includes subdirectories)
userVar_RegionResolution  <- "broadRegion"
outdir                    <- "cts_controlInteractomes/Results/"

####################################################################################### MAIN ####

### 1. Generate master PPIN from HuRI

## 1a. Import HuRI PPIs from file
HuRI <- read.delim(file = userFile_HuRI,                                                                    # Read HuRI.tsv from file
                   header = TRUE,
                   sep = "\t",
                   quote = "",
                   colClasses = "character")
HuRI <- HuRI[,1:2]                                                                                          # Remove junk columns

## 1b. In HuRI PPI list, convert Ensembl to gene symbols
# Get mappings from biomaRt
BMres <- fetchBM(unlist(HuRI)[!duplicated(unlist(HuRI))], file = userFile_biomaRt_HuRI)
# Convert first interaction partners
HuRI <- merge(x = HuRI,                                                                                     # Convert first column of HuRI using biomaRt
              y = BMres,
              by.x = "Ensembl_gene_id_a",
              by.y = "ensembl_gene_id",
              all = TRUE)
# Convert second interaction partners
HuRI <- merge(x = HuRI,                                                                                     # Convert second column of HuRI using biomaRt
              y = BMres,
              by.x = "Ensembl_gene_id_b",
              by.y = "ensembl_gene_id",
              all = FALSE)

## 1c. Remove self-interactions
HuRI <- HuRI[-(which(HuRI$hgnc_symbol.x == HuRI$hgnc_symbol.y)), ]                                          # Remove rows in which both elements are the same


### 2. Create PPINs by picking a subset of HuRI PPIs:
###    PPIs to be selected are those corresponding to pairs of genes which are expressed

## 2a. Import TPMs from file (GTEx data)
GTExTPMs <- read.delim(file = userFile_GTExTPMs,                                                            # Read GTEx TPMs from file
                       header = TRUE,
                       skip = 2,
                       sep = "\t",
                       quote = "",
                       check.names = FALSE)
# Accept only desired samples
GTExTPMs <- GTExTPMs[,c(1,which(grepl("Brain", names(GTExTPMs))))]                                          # Accept only column 1 and columns in which the column name contains "Brain" - change as desired

## 2b. Convert region labels to broad regions
# Store gene IDs separately from TPMs for now
temp <- GTExTPMs[1]                                                                                         # Store Ensembl IDs in a temporary variable
GTExTPMs <- GTExTPMs[-1]                                                                                    # Remove Ensembl IDs
# Read region mapping dictionary from file (needs to be manually curated)
RegionDict <- read.delim(file = userFile_RegionMapping,                                                     # Read the region mapping file (Column 1: region labels as in data; Column 2: precise region; Column 3: broad region)
                         header = TRUE,
                         sep = "\t",
                         quote = "")
RegionDict <- setNames(RegionDict[[userVar_RegionResolution]], RegionDict$label)                            # Create the region mapping dictionary from the region mapping file - choose between $broadRegion or $preciseRegion
# Convert labels to broad regions
colnames(GTExTPMs) <- as.character(RegionDict[colnames(GTExTPMs)])                                          # Convert region labels to broad region names
# Aggregate columns with the same name (samples from the same broad region)
GTExTPMs <- as.data.frame(do.call(cbind,                                                                    # Collapse columns corresponding to the same broad region
                                  by(t(GTExTPMs),                                                           # by transposing,
                                     INDICES = names(GTExTPMs),                                             # subsetting by broad regions,
                                     FUN = colMeans)))                                                      # and taking means across region labels
# Restore gene IDs column
GTExTPMs <- cbind(temp,                                                                                     # Add the Ensembl ID column back into the data
                  GTExTPMs)
# Remove Ensembl version numbers from Ensembl IDs
GTExTPMs[,1] <- sub(pattern = "(.+?)\\..+",                                                                 # Look for the pattern: at least one of any character (lazy), then a dot, then at least one of any character
                    replacement = "\\1",                                                                    # Replace the pattern with everything matched up to the dot
                    GTExTPMs[,1])

# 2c. Get mappings from biomaRt
BMres <- fetchBM(GTExTPMs[,1], file = userFile_biomaRt_GTExTPMs)                                            # Query biomaRt
GTExTPMs <- merge(x = GTExTPMs,                                                                             # Convert Ensembl IDs to gene symbols using biomaRt results
                  y = BMres,
                  by.x = "Name",
                  by.y = "ensembl_gene_id",
                  all = FALSE)
GTExTPMs <- cbind(GTExTPMs[1],                                                                              # Rearrange columns into the order: Ensembl ID column, then
                  GTExTPMs[ncol(GTExTPMs)],                                                                 # gene symbol column, then
                  GTExTPMs[2:(ncol(GTExTPMs)-1)])                                                           # TPMs

## 2d. Accept only expressed genes
# Mark genes with TPM >= 1 for subsetting
GTExExpressed <- cbind(GTExTPMs[1:2],                                                                       # Record Ensembl ID, gene symbol,
                       as.data.frame(GTExTPMs[,-(1:2)] >= 1))                                               # and whether TPMs are equal to or greater than 1 or not
# Restore column 1:2 class: character
GTExExpressed[,1] <- as.character(GTExExpressed[,1])
GTExExpressed[,2] <- as.character(GTExExpressed[,2])
# Remove genes with TPM < 1 (unmarked for subsetting) in all samples
GTExExpressed <- GTExExpressed[which(rowSums(GTExExpressed[,3:ncol(GTExExpressed)]) != 0), ]                # Remove rows which sum to 0 (meaning all elements are FALSE)


# 2eii. Parallel option (fast) without progress bar
GTExExpressedList <-
    foreach (r = colnames(GTExExpressed[-(1:2)])) %dopar% {                                                 # Iterate through broad regions (all column names except for Ensembl ID and gene symbol columns)
        foreach (i = 1:nrow(GTExExpressed), .combine = 'c') %do% {                                          # Iterate through genes
            if (GTExExpressed[[r]][i] == TRUE) {                                                            # If the gene is expressed in the broad region (element is TRUE)
                GTExExpressed[i,1]                                                                          # then add the Ensembl ID of that gene to a new list of expressed genes in the broad region
            }
        }
    }
names(GTExExpressedList) <- colnames(GTExExpressed[-(1:2)])                                                 # Name each list of expressed genes with the name of the broad region

## 2f. Deduce PPIs from expressed genes in each region

# 2fii. Parallel option (faster)
GTExPPINs <-
    foreach (r = names(GTExExpressedList)) %dopar% {                                                        # Iterate through broad regions, and make a new list of PPIs for broad region containing HuRI interaction pairs
        HuRI[which(HuRI$Ensembl_gene_id_a %in% GTExExpressedList[[r]] &                                     # where the first interaction partner is expressed in the broad region
                       HuRI$Ensembl_gene_id_b %in% GTExExpressedList[[r]]),                                 # and the second interaction partner is expressed in the broad region
             3:4] %>%                                                                                       # recording just the pairs of gene symbols
            filter(hgnc_symbol.x != "") %>%                                                                 # remove interactions with only one interacting partner
            filter(hgnc_symbol.y != "") %>%
            unique()
    }
names(GTExPPINs) <- names(GTExExpressedList)                                                                # Name each list of PPIs with the name of the broad regions


### 3. Create TRNs

## 3a. Read and parse FANTOM5 GGI files
# Get the list of FANTOM5 GGI files
FANTOM5Regions <- list.files(path = userDir_FANTOM5,                                                        # Find all files within the specified directory
                             pattern = ".txt",                                                              # with the ".txt" pattern in the filename
                             recursive = TRUE)                                                              # and including subdirectories
# Trim file names and save region labels
FANTOM5Regions <- data.frame(file = FANTOM5Regions,                                                         # Record filepaths
                             label = sub(pattern = ".+\\/(.+)\\.txt",                                       # and region labels (everything in filenames before the ".txt")
                                         replacement = "\\1",
                                         FANTOM5Regions),
                             stringsAsFactors = FALSE)
# Read files in the list to generate GGIs


# 3aii. Parallel option (fast) without progress bar
FANTOM5GGIs <-
    foreach (i = 1:length(FANTOM5Regions$file)) %dopar% {                                                   # Iterate through FANTOM5 files
        read.delim(file = paste(userDir_FANTOM5,                                                            # Read each file and save the GGIs
                                FANTOM5Regions$file[i],
                                sep = ""),
                   header = FALSE,
                   sep = "\t",
                   quote = "",
                   colClasses = "character")
    }
names(FANTOM5GGIs) <- FANTOM5Regions$label                                                                  # Restore the region labels to the GGIs
for (i in 1:length(FANTOM5Regions$file)) {                                                                  # Iterate through FANTOM5 files again
    colnames(FANTOM5GGIs[[FANTOM5Regions$label[i]]]) <- c("hgnc_symbol.x",                                  # Assign column names: first interaction partner (gene symbol),
                                                          "hgnc_symbol.y",                                  # second interaction partner (gene symbol)
                                                          "edge_weight")                                    # edge weight value
    # Remove edge weight column
    FANTOM5GGIs[[FANTOM5Regions$label[i]]] <- FANTOM5GGIs[[FANTOM5Regions$label[i]]][,-3]                 # Remove the edge weight value columns
}


# 3bii. Parallel option (fast) without progress bar
FANTOM5GGIs <-
    foreach (label = names(FANTOM5GGIs)) %dopar% {                                                               # Iterate through GGIs
        FANTOM5GGIs[[label]][which(FANTOM5GGIs[[label]]$hgnc_symbol.x %in% GTExExpressed$hgnc_symbol &           # Accept GGIs in which the first interaction partner is expressed in the broad region
                                       FANTOM5GGIs[[label]]$hgnc_symbol.y %in% GTExExpressed$hgnc_symbol),] %>%  # and the second interaction partner is also expressed in the broad region
            filter(hgnc_symbol.x != "") %>%                                                                      # remove interactions with only one interacting partner
            filter(hgnc_symbol.y != "") %>%
            unique()
    }
names(FANTOM5GGIs) <- FANTOM5Regions$label                                                                  # Restore the region labels to the GGIs

## 3c. Aggregate GGIs into broad regions using the region dictionary
FANTOM5TRNs <- list()
progbar <- pbinit(length(names(FANTOM5GGIs)))
for (label in names(FANTOM5GGIs)) {                                                                                # Iterate through GGIs
    FANTOM5TRNs[[as.character(RegionDict[label])]] <- unique(rbind(FANTOM5TRNs[[as.character(RegionDict[label])]], # Make a list of GGIs for each broad region, removing duplicates
                                                                   FANTOM5GGIs[[label]]))
    setTxtProgressBar(progbar, which(names(FANTOM5GGIs) == label))
}


### 4. Assemble Inxs

## 4a. Get a list of all regions covered by at least one PPIN or TRN
InxsList <- unique(c(names(GTExPPINs),                                                                      # Generate a list of broad regions for which we made a PPIN
                     names(FANTOM5TRNs)))                                                                   # and/or a TRN

## 4b. Generate Inxs from PPINs and TRNs describing the same region

# 4bii. Parallel option (fast) without progress bar
Inxs <-
    foreach (r = InxsList) %dopar% {                                                                        # Iterate through broad regions
        unique(rbind(GTExPPINs[[r]],                                                                        # Make an Inx by combining PPIs
                     FANTOM5TRNs[[r]]))                                                                     # and GGIs, removing duplicates
    }
names(Inxs) <- InxsList                                                                                     # Restore the broad region names to the Inxs






## 5i. Convert Inxs' gene symbols back to ensembl ids, save again
BMres <- rbind(read.table(file = userFile_biomaRt_GTExTPMs),
               read.table(file = userFile_biomaRt_HuRI))
Inxs2 <- foreach(r = names(Inxs),                                                                           # Convert to ensembl ids, keep ensembl and symbols
                 .final = function(r) setNames(r, names(Inxs))) %dopar% {
                     Inxs[[r]] %>%
                         inner_join(y = BMres,
                                    by = c("hgnc_symbol.x" = "hgnc_symbol")) %>%
                         inner_join(y = BMres,
                                    by = c("hgnc_symbol.y" = "hgnc_symbol")) %>%
                         unique()
                 }

FANTOM5TPMs <- read.delim(file = userFile_FANTOM5TPMs,                                                      # Read FANTOM5 TPMs from file
                          header = TRUE,
                          sep = "\t",
                          quote = "",
                          #colClasses = "character",
                          skip = 1837)                                                                      # Skip junk rows
names(FANTOM5TPMs) <- names(FANTOM5TPMs) %>%                                                                # Clean variable names
    sub(pattern = "tpm\\.",
        replacement = "") %>%
    sub(pattern = "donor.+",
        replacement = "") %>%
    sub(pattern = "\\.CNhs.+",
        replacement = "") %>%
    gsub(pattern = "\\.20\\.\\.20",
         replacement = "_-_") %>%
    gsub(pattern = "\\.20",
         replacement = "_") %>%
    gsub(pattern = "\\.2c",
         replacement = "") %>%
    gsub(pattern = "\\.",
         replacement = "") %>%
    gsub(pattern = "[ ._-]+$",
         replacement = "") %>%
    tolower()

temp <- FANTOM5TPMs[3]                                                                                      # Store observation names for later
FANTOM5TPMs <- FANTOM5TPMs[-(1:7)]                                                                          # Remove observation metadata
#FANTOM5TPMs <- FANTOM5TPMs[which(names(FANTOM5TPMs) %in% names(RegionDict))]       # old
#FANTOM5TPMs <- FANTOM5TPMs[sapply(X = names(RegionDict), FUN = function(X) grep(x = names(FANTOM5TPMs), pattern = X)) %>% unlist() %>% as.numeric()] # old
FANTOM5TPMs <- foreach (r = names(RegionDict),                                                                      # Loop to accept only variables which are in the region dictionary
                .final = function(r) setNames(r, as.character(RegionDict)
                )
) %do% {
    grep(x = names(FANTOM5TPMs), pattern = r) %>%
        unlist() %>%
        as.numeric() %>%
        dplyr::select(.data = FANTOM5TPMs)
}

FANTOM5TPMs <- foreach (r = unique(names(FANTOM5TPMs)),                                                                     # Loop to aggregate variables corresponsing to the same region
                .final = function(r) setNames(r, unique(names(FANTOM5TPMs)
                )
                )
) %do% {
    do.call(what = cbind,
            args = FANTOM5TPMs[which(names(FANTOM5TPMs) == r)]) %>%
        rowMeans()
}
#colnames(FANTOM5TPMs) <- as.character(RegionDict[colnames(FANTOM5TPMs)])                                    # Convert variable names by the region dictionary
#FANTOM5TPMs <- as.data.frame(sapply(X = FANTOM5TPMs,                                                        # Convert the table cell classes to numeric
#                                    FUN = as.numeric))
#FANTOM5TPMs <- as.data.frame(do.call(cbind,                                                                 # Collapse columns corresponding to the same broad region
#                                     by(t(FANTOM5TPMs),                                                     # by transposing,
#                                        INDICES = names(FANTOM5TPMs),                                       # subsetting by broad regions,
#                                        FUN = colMeans)))
FANTOM5TPMs <- cbind(temp, FANTOM5TPMs)                                                                     # Restore observation names
FANTOM5TPMs <- FANTOM5TPMs %>%
    filter(grepl(pattern = "ENST",                                                                          # Accept only observations for which ensembl transcript id is known
                 x = description)) %>%
    mutate(description = regmatches(x = description,                                                        # Clean accepted observation names
                                    m = regexpr(text = description,
                                                pattern = "ENST\\d+")))
FANTOM5TPMs <- FANTOM5TPMs[(which(rowSums(FANTOM5TPMs[2:ncol(FANTOM5TPMs)],
                                          na.rm = T) != 0)),]                        # Reject observations for which 0 TPMs were recorded in all regions
BMres <- getBM(attributes = c("ensembl_gene_id",                                                            # Use biomaRt to convert ENST to ENSG
                              "ensembl_transcript_id"),
               filters = "ensembl_transcript_id",
               values = FANTOM5TPMs$description,
               mart = useMart("ensembl",
                              dataset = "hsapiens_gene_ensembl",
                              host = "www.ensembl.org"))
FANTOM5TPMs <- FANTOM5TPMs %>%                                                                              # Join ENSG ids to the table
    inner_join(y = BMres,
               by = c("description" = "ensembl_transcript_id"))



Inxs3 <- foreach(r = names(Inxs),                                                                           # Keep only ensembl ids, collapse list
                 .final = function(r) setNames(r, names(Inxs))) %dopar% {
                     l <- Inxs2[[r]] %>%
                         dplyr::select(ensembl_gene_id.x,
                                       ensembl_gene_id.y) %>%
                         unlist() %>%
                         unique() %>%
                         as.character() %>%
                         as_tibble()
                     try({                                                                                  # Try to get GTEx TPM data
                         l <- l %>%
                             left_join(y = GTExTPMs %>% dplyr::select(Name, r),
                                       by = c("value" = "Name") )%>%
                             filter(!duplicated(value))
                     }, silent = T)
                     try({                                                                                  # Try to get FANTOM5 TPM data
                         l <- l %>%
                             left_join(y = FANTOM5TPMs %>% dplyr::select(ensembl_gene_id, r),
                                       by = c("value" = "ensembl_gene_id")) %>%
                             filter(!duplicated(value))
                     }, silent = T)
                     if (ncol(l) > 2) {
                         l <- l %>%                                                                             # Scale (within data source) and then average (if more than one source)
                         transmute(value = value,
                                   expr = removeBatchEffect(x = l[2:ncol(l)], batch = LETTERS[1:ncol(l)-1]) %>% rowMeans(na.rm = T))
                     }
                     l
                 }

dir.create(path = outdir, showWarnings = FALSE)
dir.create(path = paste0(outdir, userVar_RegionResolution), showWarnings = FALSE)


for (r in names(Inxs)) {
    write.table(x = Inxs3[[r]],                                                                             # Save all interactions in Inxs for each broad region to file
                file = paste0(outdir,
                             userVar_RegionResolution,
                             "/",
                             as.character(r),
                             ".tsv"),
                sep = "\t")
}


## 5i. Visualise the Inxs (extremely slow)
# say("Visualising the Inxs...")
# for (r in 1:length(names(Inxs))) {
#     png(paste("Results/Inxs/",
#               names(Inxs)[r],
#               ".png",
#               sep = ""),
#         width = 1000,
#         height = 1000)
#     g <- graph_from_data_frame(d = Inxs[[r]])
#     g <- ggraph(g, layout = "kk") +
#         geom_edge_link(alpha = 0.3) +
#         geom_node_point(alpha = 1) +
#         # geom_node_text(aes(label = name), repel = TRUE) +
#         ggtitle(paste("Inx for",
#                       names(Inxs)[r]))
#     print(g)
#     dev.off()
# }

### 6. CLEAN UP WORKSPACE

## 7b. Clean up parallelisation environment
stopImplicitCluster()


################################################################################### METADATA ####
#
# zeb_cts.R
# Version: 1.0.1 (2020-02-27)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Quantification and normalisation of counts for zebrafish RNA-seq.
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
##################################################################################### INPUTS ####
#
# [1] zeb/Data/kallisto/
#     Directory containing kallisto quant results.
#
#################################################################################### OUTPUTS ####
#
# [1] zeb/Data/processedCts/
#     Normalised counts
#
# [2] zeb/Data/finalCts/finalCts.tsv
#     Cleaned and combined counts ready for submission to summarizedExperiment.
#
# [3] zeb/Data/finalCts/finalCts.se.Rdata
#     Cleaned and combined counts and adjacent metadata ready for DESeq2 analysis.
#
# [4] zeb/Data/finalCts/colData.tsv
#     Metadata.
#
########################################################################## FILE REQUIREMENTS ####

require(data.table)
require(foreach)
require(doParallel)
require(dplyr)
require(tidyr)
require(SummarizedExperiment)

kallistoDir <- "zeb/Data/kallisto/"
processedCtsDir <- "zeb/Data/processedCts/"
finalColData <- "zeb/Data/finalCts/colData.tsv"
finalCts <- "zeb/Data/finalCts/finalCts.tsv"
finalSE <- "zeb/Data/finalCts/finalCts.se.Rdata"

registerDoParallel(detectCores() - 2)

####################################################################################### MAIN ####

counts <- foreach (dir = list.files(kallistoDir),
		   .final = function(dir) setNames(dir, list.files(kallistoDir))) %dopar% {
		if (grepl(".fastq.gz-counts", dir)) {
			fread(paste(kallistoDir,
                                    dir,
				    "/abundance.tsv",
                	            sep = "")) %>%
				transmute(target_id = target_id,
					  counts = round(est_counts * eff_length / length),
                	                  id = paste(strsplit(dir, "_")[[1]][1],
                       	                             strsplit(dir, "_")[[1]][2]))
		}
}
foreach (file = names(counts)) %dopar% {
        if (grepl(".fastq.gz-counts", file)) {
		fwrite(x = counts[[file]],
       		       file = paste(processedCtsDir,
                                    file,
                            	    ".csv",
                                    sep = ""))
	}
}

counts2 <- tibble(target_id = "",
                 counts = "",
                 specimen = "")

for (file in list.files(processedCtsDir)) {
	counts2 <- counts2 %>%
		    rbind(fread(paste(processedCtsDir,
                                          file,
                                          sep = "")),
                          use.names = F)
}

counts2 <- counts2 %>%
	    filter(target_id != "") %>%
	    spread(specimen, counts)
counts2 <- counts2 %>%
            mutate(target_id = sub(x = counts2$target_id, pattern = "(.+)\\..+", replacement = "\\1"))

fwrite(x = counts2, file = paste0(finalCts), sep = "\t")


### 1. Read and wrangle raw counts file
counts3 <- list()
colData <- list()

## 1a. Wrangle Full model - counts
counts3$Full <- fread(paste0(finalCts))                                             # Import full merged counts data
counts3$Full <- as.data.frame(counts3$Full)
rownames(counts3$Full) <- counts3$Full[,1]                                        # Move Geneid into row names
counts3$Full <- counts3$Full[2:length(colnames(counts3$Full))]

## 1b. Wrangle Full model - colData
colData$Full <- data.frame(filename = colnames(counts3$Full),                    # Make full colData, ignore first column name "Geneid"
                          stringsAsFactors = F)
colData$Full$filename <- regmatches(x = colData$Full$filename,
                                   regexpr(pattern = "[[:upper:]].+?\\d{3}",    # Clean filename names
                                           text = colData$Full$filename))
colData$Full <- colData$Full %>%
    separate(col = filename,                                                    # Separate the filename name column into part and specimen
             into = c("part", "specimen")) %>%
    mutate(tert1 = case_when(specimen == 234 | specimen == 253 | specimen == 255 | specimen == 257 | specimen == 270 ~ "wt",
                             specimen == 235 | specimen == 251 | specimen == 258 | specimen == 264 | specimen == 271 ~ "het",
                             specimen == 254 | specimen == 256 | specimen == 273 ~ "dko")) %>% # Populate tert1 data
    mutate(sample = paste("NH_FLI", specimen, part, tert1, sep = "_")) %>%      # Clean sample names
    transmute(sample = sample,                                                  # Rearrange fields
              part = part,
              specimen = specimen,
              sex = "male",                                                     # Populate common fields
              ageMonths = 3,
              tert1 = tert1)

colnames(counts3$Full) <- colData$Full$sample                                    # Apply cleaned sample names back to counts matrix

### 2. Wrangle subsetted models
## 2a. Subset Full model by part
for (partSel in c("Brain", "Liver", "Muscle", "Skin")) {                        # By part
    counts3[[partSel]] <- counts3$Full %>%
        dplyr::select(contains(partSel))                                               # subset counts
    colData[[partSel]] <- colData$Full %>%
        filter(part == partSel)                                                 # and colData
}

## 2b. Subset Full model by genotype
for (tert1Sel in c("wt", "het", "dko")) {                                       # By genotype
    counts3[[tert1Sel]] <- counts3$Full %>%
        dplyr::select(contains(tert1Sel))                                              # subset counts
    colData[[tert1Sel]] <- colData$Full %>%
        filter(tert1 == tert1Sel)                                               # and colData
}

### 3. Save to file
if (identical(names(counts3), names(colData))) {                                 # Check to see if the subsetting worked properly
    se <- list()                                                                # Initialise a list of ses
    for (name in names(counts3)) {
        se[[name]] <- SummarizedExperiment(assays = as.matrix(counts3[[name]]),  # Create an se for each subset
                                           colData = colData[[name]])
    }
    save(se, file = finalSE)                                                 # Save summarizedExperiment list to file
}

colData_out <- colData$Full %>%
    mutate(sample = sub(x = colData$Full$sample, pattern = "NH_FLI_(\\d{3})_(.+)_.+", replacement = "\\2 \\1"))

fwrite(colData_out, finalColData, sep = "\t")

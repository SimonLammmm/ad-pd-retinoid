################################################################################### METADATA ####
#
# cts_ROSMAP_blood_RNAseq_kallisto.sh
# Version: 1.0.0 (2020-08-11)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Quantification of counts for ROSMAP blood RNA-seq.
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
# [1] cts_ROSMAP_blood_RNAseq/Data/Sample_XXX_ENDY.fastq.gz-trimmed.fastq.gz
#     Trimmed fastq reads.
#
##################################################################################### INPUTS ####
#
# [1] kallisto/homo_sapiens/transcriptome.idx
#     Homo sapiens transcriptome index for kallisto.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_blood_RNAseq/Data/Sample_XXX_END1.fastq.gz-trimmed.fastq.gz-kallisto_out/
#     Directory containing kallisto quant results
#
####################################################################################### MAIN ####
for FILE1 in cts_ROSMAP_blood_RNAseq/Data/*END1.fastq.gz-trimmed.fastq.gz
do
	FILE2=`echo $FILE1 | sed 's/END1/END2/g'`
	mkdir $FILE1-kallisto_out
	kallisto quant \
	-i kallisto/homo_sapiens/transcriptome.idx \
	-o $FILE1-kallisto_out \
	$FILE1 $FILE2
done

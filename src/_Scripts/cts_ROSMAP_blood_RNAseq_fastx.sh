################################################################################### METADATA ####
#
# cts_ROSMAP_blood_RNAseq_fastx.sh
# Version: 1.0.0 (2020-07-27)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Trimming of ROSMAP blood RNA-seq raw reads with FASTX.
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
########################################################################## FILE REQUIREMENTS ####
#
# [1] cts_ROSMAP_blood_RNAseq/Data/Sample_XXX_ENDY.fastq.gz
#     https://www.synapse.org/#!Synapse:syn22024498
#     Approval required.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_blood_RNAseq/Data/Sample_XXX_ENDY.fastq.gz-trimmed.fastq.gz
#     Trimmed fastq reads.
#
####################################################################################### MAIN ####
for FILE in cts_ROSMAP_blood_RNAseq/Data/*.fastq.gz
do
	echo "processing $FILE"
	gunzip -c $FILE | fastx_trimmer -f 16 -z -o "$FILE-trimmed.fastq.gz"
	echo "\n"
done

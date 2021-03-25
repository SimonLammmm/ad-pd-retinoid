################################################################################### METADATA ####
#
# cts_ROSMAP_microglia_single_cell_RNA_seq_rename.sh
# Version: 1.0.0 (2020-11-12)
#
# Author: Simon Lam
# Institution: King's College London
# Contact: simon.1.lam@kcl.ac.uk
#
# Description:
# Quantification of counts for ROSMAP microglia single-cell RNA-seq.
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
# [1] cts_ROSMAP_microglia_single_cell_RNA_seq/Data/Microglia_MO_XXXYY_LZZZ_AB_001.fastq.gz
#     https://www.synapse.org/#!Synapse:syn21438358
#     Cell Ranger fastq raw reads.  Approval required.  Follow the  renaming instructions  on the
#     webpage or execute this script.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_microglia_single_cell_RNA_seq/Data/Microglia_MO_XXXYY_S1_LZZZ_AB_001.fastq.gz
#     https://www.synapse.org/#!Synapse:syn21438358
#     Cell Ranger fastq raw reads, correctly renamed for Cell Ranger count.
#
####################################################################################### MAIN ####
CR_IN='cts_ROSMAP_microglia_single_cell_RNA_seq/Data'
for FILE in $CR_IN/*
do
	NEWNAME=`echo $FILE | sed s/L00/S1_L00/g`
	mv $FILE $NEWNAME
done

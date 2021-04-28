################################################################################### METADATA ####
#
# cts_ROSMAP_microglia_single_cell_RNA_seq_cellranger.sh
# Version: 1.0.0 (2020-08-19)
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
#     Approval required.
#
# [2] cellranger/refdata-gex-GRCh38-2020-A/
#     https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
#     Sign-up required at:
#     https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_microglia_single_cell_RNA_seq/Results/CellRanger/
#     Directories containing Cell Ranger counts
#
####################################################################################### MAIN ####
OLD_PWD=`pwd`
mkdir cts_ROSMAP_microglia_single_cell_RNA_seq/Results/CellRanger/
cd cts_ROSMAP_microglia_single_cell_RNA_seq/Results/CellRanger/
for SAMPLE in Microglia_MO_AD1 Microglia_MO_AD2 Microglia_MO_AD3 Microglia_MO_AD4 Microglia_MO_AD51 Microglia_MO_AD6 Microglia_MO_AD71 Microglia_MO_AD8 Microglia_MO_AD9 Microglia_MO_MCI1 Microglia_MO_MCI2 Microglia_MO_MCI3 Microglia_MO_MCI4
do
	cellranger count --id=$SAMPLE-counts --fastqs=../../Data --transcriptome=../../../cellranger/refdata-gex-GRCh38-2020-A --sample=$SAMPLE
done
cd $OLD_PWD

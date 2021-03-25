################################################################################### METADATA ####
#
# cts_ROSMAP_single_nucleus_RNA_seq_rename.sh
# Version: 1.0.0 (2020-10-29)
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
# [1] cts_ROSMAP_single_nucleus_RNA_seq/Data/MFC-B1-SX-CdxY-pADZ_CD.fastq.gz
#     https://www.synapse.org/#!Synapse:syn16780177
#     Cell Ranger fastq raw reads.  Approval required.  Follow the  renaming instructions  on the
#     webpage or execute this script.
#
#################################################################################### OUTPUTS ####
#
# [1] cts_ROSMAP_single_nucleus_RNA_seq/Data/MFC-B1-SX-CdxY-pADZ_SA_LBBB_CD_001.fastq.gz
#     https://www.synapse.org/#!Synapse:syn16780177
#     Cell Ranger fastq raw reads, correctly renamed for Cell Ranger count.
#
####################################################################################### MAIN ####
CR_INDEX_IN='cts_ROSMAP_single_nucleus_RNA_seq/Data/Index'
CR_FASTQ_IN='cts_ROSMAP_single_nucleus_RNA_seq/Data/fastq'
CR_BATCH2_IN='cts_ROSMAP_single_nucleus_RNA_seq/Data/fastq/batch2'
CR_OUT='cts_ROSMAP_single_nucleus_RNA_seq/Data'
mv $CR_INDEX_IN/MFC-B1-S1-Cdx1-pAD0_L006_I1_001.fastq.gz $CR_OUT/MFC-B1-S1-Cdx1-pAD0_S1_L006_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S1-Cdx1-pAD0_L006_R1_001.fastq.gz $CR_OUT/MFC-B1-S1-Cdx1-pAD0_S1_L006_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S1-Cdx1-pAD0_L006_R2_001.fastq.gz $CR_OUT/MFC-B1-S1-Cdx1-pAD0_S1_L006_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S2-Cdx1-pAD0_L007_I1_001.fastq.gz $CR_OUT/MFC-B1-S2-Cdx1-pAD0_S4_L007_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S2-Cdx1-pAD0_L007_R1_001.fastq.gz $CR_OUT/MFC-B1-S2-Cdx1-pAD0_S4_L007_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S2-Cdx1-pAD0_L007_R2_001.fastq.gz $CR_OUT/MFC-B1-S2-Cdx1-pAD0_S4_L007_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S3-Cdx1-pAD1_L007_I1_001.fastq.gz $CR_OUT/MFC-B1-S3-Cdx1-pAD1_S1_L007_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S3-Cdx1-pAD1_L007_R1_001.fastq.gz $CR_OUT/MFC-B1-S3-Cdx1-pAD1_S1_L007_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S3-Cdx1-pAD1_L007_R2_001.fastq.gz $CR_OUT/MFC-B1-S3-Cdx1-pAD1_S1_L007_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S4-Cdx1-pAD1_L007_I1_001.fastq.gz $CR_OUT/MFC-B1-S4-Cdx1-pAD1_S5_L007_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S4-Cdx1-pAD1_L007_R1_001.fastq.gz $CR_OUT/MFC-B1-S4-Cdx1-pAD1_S5_L007_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S4-Cdx1-pAD1_L007_R2_001.fastq.gz $CR_OUT/MFC-B1-S4-Cdx1-pAD1_S5_L007_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S5-Cdx4-pAD0_L006_I1_001.fastq.gz $CR_OUT/MFC-B1-S5-Cdx4-pAD0_S2_L006_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S5-Cdx4-pAD0_L006_R1_001.fastq.gz $CR_OUT/MFC-B1-S5-Cdx4-pAD0_S2_L006_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S5-Cdx4-pAD0_L006_R2_001.fastq.gz $CR_OUT/MFC-B1-S5-Cdx4-pAD0_S2_L006_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S6-Cdx4-pAD0_L007_I1_001.fastq.gz $CR_OUT/MFC-B1-S6-Cdx4-pAD0_S6_L007_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S6-Cdx4-pAD0_L007_R1_001.fastq.gz $CR_OUT/MFC-B1-S6-Cdx4-pAD0_S6_L007_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S6-Cdx4-pAD0_L007_R2_001.fastq.gz $CR_OUT/MFC-B1-S6-Cdx4-pAD0_S6_L007_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S7-Cdx4-pAD1_L006_I1_001.fastq.gz $CR_OUT/MFC-B1-S7-Cdx4-pAD1_S3_L006_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S7-Cdx4-pAD1_L006_R1_001.fastq.gz $CR_OUT/MFC-B1-S7-Cdx4-pAD1_S3_L006_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S7-Cdx4-pAD1_L006_R2_001.fastq.gz $CR_OUT/MFC-B1-S7-Cdx4-pAD1_S3_L006_R2_001.fastq.gz
mv $CR_INDEX_IN/MFC-B1-S8-Cdx4-pAD1_L007_I1_001.fastq.gz $CR_OUT/MFC-B1-S8-Cdx4-pAD1_S7_L007_I1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S8-Cdx4-pAD1_L007_R1_001.fastq.gz $CR_OUT/MFC-B1-S8-Cdx4-pAD1_S7_L007_R1_001.fastq.gz
mv $CR_FASTQ_IN/MFC-B1-S8-Cdx4-pAD1_L007_R2_001.fastq.gz $CR_OUT/MFC-B1-S8-Cdx4-pAD1_S7_L007_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-10-Cog1-Path0_L005_I1_001.fastq.gz $CR_OUT/MFC-B2-10-Cog1-Path0_S1_L005_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-10-Cog1-Path0_L005_R1_001.fastq.gz $CR_OUT/MFC-B2-10-Cog1-Path0_S1_L005_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-10-Cog1-Path0_L005_R2_001.fastq.gz $CR_OUT/MFC-B2-10-Cog1-Path0_S1_L005_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-11-Cog1-Path1_L004_I1_001.fastq.gz $CR_OUT/MFC-B2-11-Cog1-Path1_S1_L004_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-11-Cog1-Path1_L004_R1_001.fastq.gz $CR_OUT/MFC-B2-11-Cog1-Path1_S1_L004_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-11-Cog1-Path1_L004_R2_001.fastq.gz $CR_OUT/MFC-B2-11-Cog1-Path1_S1_L004_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-12-Cog1-Path1_L005_I1_001.fastq.gz $CR_OUT/MFC-B2-12-Cog1-Path1_S1_L005_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-12-Cog1-Path1_L005_R1_001.fastq.gz $CR_OUT/MFC-B2-12-Cog1-Path1_S1_L005_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-12-Cog1-Path1_L005_R2_001.fastq.gz $CR_OUT/MFC-B2-12-Cog1-Path1_S1_L005_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-13-Cog4-Path0_L004_I1_001.fastq.gz $CR_OUT/MFC-B2-13-Cog4-Path0_S1_L004_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-13-Cog4-Path0_L004_R1_001.fastq.gz $CR_OUT/MFC-B2-13-Cog4-Path0_S1_L004_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-13-Cog4-Path0_L004_R2_001.fastq.gz $CR_OUT/MFC-B2-13-Cog4-Path0_S1_L004_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-14-Cog4-Path0_L005_I1_001.fastq.gz $CR_OUT/MFC-B2-14-Cog4-Path0_S1_L005_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-14-Cog4-Path0_L005_R1_001.fastq.gz $CR_OUT/MFC-B2-14-Cog4-Path0_S1_L005_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-14-Cog4-Path0_L005_R2_001.fastq.gz $CR_OUT/MFC-B2-14-Cog4-Path0_S1_L005_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-15-Cog4-Path1_L004_I1_001.fastq.gz $CR_OUT/MFC-B2-15-Cog4-Path1_S1_L004_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-15-Cog4-Path1_L004_R1_001.fastq.gz $CR_OUT/MFC-B2-15-Cog4-Path1_S1_L004_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-15-Cog4-Path1_L004_R2_001.fastq.gz $CR_OUT/MFC-B2-15-Cog4-Path1_S1_L004_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-16-Cog4-Path1_L005_I1_001.fastq.gz $CR_OUT/MFC-B2-16-Cog4-Path1_S1_L005_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-16-Cog4-Path1_L005_R1_001.fastq.gz $CR_OUT/MFC-B2-16-Cog4-Path1_S1_L005_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-16-Cog4-Path1_L005_R2_001.fastq.gz $CR_OUT/MFC-B2-16-Cog4-Path1_S1_L005_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-9-Cog1-Path0_L004_I1_001.fastq.gz $CR_OUT/MFC-B2-9-Cog1-Path0_S1_L004_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-9-Cog1-Path0_L004_R1_001.fastq.gz $CR_OUT/MFC-B2-9-Cog1-Path0_S1_L004_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B2-9-Cog1-Path0_L004_R2_001.fastq.gz $CR_OUT/MFC-B2-9-Cog1-Path0_S1_L004_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-17-Cog1-Path0_L001_I1_001.fastq.gz $CR_OUT/MFC-B3-17-Cog1-Path0_S1_L001_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-17-Cog1-Path0_L001_R1_001.fastq.gz $CR_OUT/MFC-B3-17-Cog1-Path0_S1_L001_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-17-Cog1-Path0_L001_R2_001.fastq.gz $CR_OUT/MFC-B3-17-Cog1-Path0_S1_L001_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-18-Cog1-Path0_L002_I1_001.fastq.gz $CR_OUT/MFC-B3-18-Cog1-Path0_S1_L002_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-18-Cog1-Path0_L002_R1_001.fastq.gz $CR_OUT/MFC-B3-18-Cog1-Path0_S1_L002_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-18-Cog1-Path0_L002_R2_001.fastq.gz $CR_OUT/MFC-B3-18-Cog1-Path0_S1_L002_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-19-Cog1-Path1_L001_I1_001.fastq.gz $CR_OUT/MFC-B3-19-Cog1-Path1_S1_L001_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-19-Cog1-Path1_L001_R1_001.fastq.gz $CR_OUT/MFC-B3-19-Cog1-Path1_S1_L001_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-19-Cog1-Path1_L001_R2_001.fastq.gz $CR_OUT/MFC-B3-19-Cog1-Path1_S1_L001_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-20-Cog1-Path1_L002_I1_001.fastq.gz $CR_OUT/MFC-B3-20-Cog1-Path1_S1_L002_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-20-Cog1-Path1_L002_R1_001.fastq.gz $CR_OUT/MFC-B3-20-Cog1-Path1_S1_L002_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-20-Cog1-Path1_L002_R2_001.fastq.gz $CR_OUT/MFC-B3-20-Cog1-Path1_S1_L002_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-21-Cog4-Path0_L001_I1_001.fastq.gz $CR_OUT/MFC-B3-21-Cog4-Path0_S1_L001_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-21-Cog4-Path0_L001_R1_001.fastq.gz $CR_OUT/MFC-B3-21-Cog4-Path0_S1_L001_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-21-Cog4-Path0_L001_R2_001.fastq.gz $CR_OUT/MFC-B3-21-Cog4-Path0_S1_L001_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-22-Cog4-Path0_L002_I1_001.fastq.gz $CR_OUT/MFC-B3-22-Cog4-Path0_S1_L002_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-22-Cog4-Path0_L002_R1_001.fastq.gz $CR_OUT/MFC-B3-22-Cog4-Path0_S1_L002_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-22-Cog4-Path0_L002_R2_001.fastq.gz $CR_OUT/MFC-B3-22-Cog4-Path0_S1_L002_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-23-Cog4-Path1_L001_I1_001.fastq.gz $CR_OUT/MFC-B3-23-Cog4-Path1_S1_L001_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-23-Cog4-Path1_L001_R1_001.fastq.gz $CR_OUT/MFC-B3-23-Cog4-Path1_S1_L001_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-23-Cog4-Path1_L001_R2_001.fastq.gz $CR_OUT/MFC-B3-23-Cog4-Path1_S1_L001_R2_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-24-Cog4-Path1_L002_I1_001.fastq.gz $CR_OUT/MFC-B3-24-Cog4-Path1_S1_L002_I1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-24-Cog4-Path1_L002_R1_001.fastq.gz $CR_OUT/MFC-B3-24-Cog4-Path1_S1_L002_R1_001.fastq.gz
mv $CR_BATCH2_IN/MFC-B3-24-Cog4-Path1_L002_R2_001.fastq.gz $CR_OUT/MFC-B3-24-Cog4-Path1_S1_L002_R2_001.fastq.gz

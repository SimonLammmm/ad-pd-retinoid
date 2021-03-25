=================================================================================== METADATA ====

README.txt
Version: 1.0.1 (2021-03-09)

Author: Simon Lam
Institution: King's College London
Contact: simon.1.lam@kcl.ac.uk

Description:
Provenance, code, and execution  instructions to ensure reproducibility  of results regarding the
human-zebrafish retinoids and sex hormones project.

====================================================================================== LEGAL ====

Copyright © 2020 King's College London

This work is licensed under the Creative Commons Attribution 4.0 International Licence. To view a
copy  of this  license,  visit http://creativecommons.org/licences/by/4.0/  or send  a letter  to
Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

Permission is  hereby granted,  free of charge, to any person  obtaining a copy  of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use,  copy, modify, merge, publish, distribute, and/or
sell copies  of the Software, and to permit persons  to whom the Software  is furnished to do so,
subject to the following conditions:

The above  copyright  notice  and this  permission notice  shall be included  with all  copies or
substantial portions of the Software.

THE SOFTWARE  IS PROVIDED  "AS IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR IMPLIED, INCLUDING
BUT NOT  LIMITED TO THE  WARRANTIES  OF MERCHANTABILITY,  FITNESS FOR A  PARTICULAR  PURPOSE, AND
NONINFRINGEMENT.  IN NO EVENT SHALL  THE AUTHORS  OR COPYRIGHT  HOLDERS BE  LIABLE FOR ANY CLAIM,
DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR  OTHERWISE, ARISING FROM,
OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF OR DEALING IN THE SOFTWARE.

================================================================== PROGRAMMATIC REQUIREMENTS ====

 1. R 3.6.2
     a. gridExtra_2.3
	 b. umap_0.2.6.0																		[1]
	 c. ggbiplot_0.55
	 d. plyr_1.8.6
	 e. writexl_1.3
	 f. piano_2.2.0																			[2]
	 g. scales_1.1.1
	 h. DESeq2_1.26.0																		[3]
	 i. SummarizedExperiment_1.16.1
	 j. DelayedArray_0.12.2
	 k. BiocParallel_1.20.1
	 l. matrixStats_0.56.0
	 m. GenomicRanges_1.38.0
	 n. GenomeInfoDb_1.22.0
	 o. IRanges_2.20.2
	 p. S4Vectors_0.24.3
	 q. hgu133acdf_2.18.0
	 r. affy_1.64.0
	 s. Biobase_2.46.0
	 t. BiocGenerics_0.32.0
	 u. R.utils_2.9.2
	 v. R.oo_1.23.0
	 w. R.methodsS3_1.8.0
	 x. tidyselect_1.1.0
	 y. tidyr_1.1.1
	 z. doParallel_1.0.15
	aa. iterators_1.0.12
	ab. ConsensusClusterPlus_1.50.0															[4]
	ac. ggplot2_3.3.2
	ad. Rtsne_0.15																			[5]
	ae. MetabolAnalyze_1.3.1																[6]
	af. gplots_3.0.4
	ag. gtools_3.8.2
	ah. ellipse_0.4.2
	ai. mvtnorm_1.1-1
	aj. mclust_5.4.6
	ak. Rmagic_2.0.3																		[7]
	al. Matrix_1.2-18
	am. tibble_3.0.3
	an. edgeR_3.28.1																		[8]
	ao. limma_3.42.0																		[9]
	ap. readxl_1.3.1
	aq. biomaRt_2.42.0																	[10,11]
	ar. foreach_1.5.0
	as. data.table_1.12.8
	at. dplyr_1.0.2

 2. Bash
    a. FASTX-Toolkit
	   http://hannonlab.cshl.edu/fastx_toolkit/download.html
    b. kallisto																			   [12]
	   https://pachterlab.github.io/kallisto/download
	c. Cell Ranger 4.0
	   https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

3. MATLAB R2020b
    a. CellFie																			   [15]
	   https://github.com/LewisLabUCSD/CellFie
	b. COBRA Toolbox																	   [13]
	   https://github.com/opencobra/cobratoolbox
	c. RAVEN Toolbox 2																	[14,15]
	   https://github.com/SysBioChalmers/RAVEN
	d. LibSBML																			   [16]
	   https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/

========================================================================== FILE REQUIREMENTS ====

 1. cellranger/refdata-gex-GRCh38-2020-A/
    https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
 2. combineCounts/Data/metadata.xlsx
    Provided.
 3. cts_controlInteractomes/FANTOM5_individual_networks/								   [17]
    http://www2.unil.ch/cbg/regulatorycircuits/FANTOM5_individual_networks.tar
 4. cts_controlInteractomes/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
    https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
 5. cts_controlInteractomes/HuRI.tsv													   [18]
    http://www.interactome-atlas.org/data/HuRI.tsv
 6. cts_controlInteractomes/RegionMapping.txt
    Provided.
 7. cts_HPA/Data/normal_tissue.tsv														   [19]
    https://www.proteinatlas.org/download/normal_tissue.tsv.zip
 8. cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc2.xlsx								   [20]
    https://ars.els-cdn.com/content/image/1-s2.0-S1064748119304142-mmc2.xlsx
 9. cts_rajkumar/Data/1-s2.0-S1064748119304142-mmc3.xlsx								   [20]
    https://ars.els-cdn.com/content/image/1-s2.0-S1064748119304142-mmc3.xlsx
10. cts_ROSMAP_blood_RNAseq/Data/														   [20]
    https://www.synapse.org/#!Synapse:syn22024498
11. cts_ROSMAP_bulk_brain_RNA_seq/Data/ROSMAP_RNAseq_FKPM_gene.tsv						   [20]
    https://www.synapse.org/#!Synapse:syn3505720
12. cts_ROSMAP_microglia_RNA_seq/Data/													   [20]
    https://www.synapse.org/#!Synapse:syn11578941
13. cts_ROSMAP_microglia_single_cell_RNA_seq/Data/										   [20]
    https://www.synapse.org/#!Synapse:syn21438358
14. cts_ROSMAP_RNA_array/Data/ROSMAP_arrayExpression_normalized.tsv						   [20]
    https://www.synapse.org/#!Synapse:syn4009614
15. cts_ROSMAP_single_nucleus_RNA_seq/Data/												   [20]
    https://www.synapse.org/#!Synapse:syn16780177
16. cts_zhengZhang/Data/GSE20295_RAW/													[21,22]
    https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE20295&format=file
17. GEM/Data/common_tasks_growth_RPMI1640.xlsx
    https://github.com/sysmedicine/phd2020/blob/master/GEM/data/Data.zip
    Data.zip/Lab/common_tasks_growth_RPMI1640.xlsx
18. iBrain2845/Data/HMRdatabase3_00.xlsx
    Provided.
19. iBrain2845/Data/iAdipocytes1850_v2.xlsx
    Provided.
20. kallisto/danio_rerio/Danio_rerio.GRCz11.cdna.all.fa									   [23]
    ftp://ftp.ensembl.org/pub/release-96/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
21. kallisto/homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa								   [23]
    ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
22. zeb/Data/fastq/
    Deposited in GEO: GSE102426, GSE102429, GSE102431, GSE102434.
23. ZebraGEM/Data/ZebraGEM2.1.xlsx														   [24]
    Provided.

============================================================================ EXECUTION ORDER ====

Processing of human brain expression counts
 1. _Scripts/cts_controlInteractomes-broad.R
 2. _Scripts/cts_controlInteractomes-precise.R
 3. _Scripts/cts_zhangZheng.R
 4. _Scripts/cts_rajkumar.R
 5. _Scripts/cts_HPA.R
 6. _Scripts/cts_ROSMAP_blood_RNAseq_fastx.sh
 7. _Scripts/kallisto_index_hsa.sh
 8. _Scripts/cts_ROSMAP_blood_RNAseq_kallisto.sh
 9. _Scripts/cts_ROSMAP_blood_RNAseq.R
10. _Scripts/cts_ROSMAP_microglia_RNA_seq_fastx.sh
11. _Scripts/cts_ROSMAP_microglia_RNA_seq_kallisto.sh
12. _Scripts/cts_ROSMAP_microglia_RNA_seq.R
13. _Scripts/cts_ROSMAP_microglia_single_cell_RNA_seq_rename.sh
14. _Scripts/cts_ROSMAP_microglia_single_cell_RNA_seq_cellranger.sh
15. _Scripts/cts_ROSMAP_microglia_single_cell_RNA_seq.R
16. _Scripts/cts_ROSMAP_single_nucleus_RNA_seq_rename.sh
17. _Scripts/cts_ROSMAP_single_nucleus_RNA_seq_cellranger.sh
18. _Scripts/cts_ROSMAP_single_nucleus_RNA_seq.R
19. _Scripts/cts_ROSMAP_bulk_brain_RNA_seq.R
20. _Scripts/cts_ROSMAP_RNA_array.R

Generation of a generalised GEM
21. _Scripts/iBrain2845.R

Combining and analysing human brain expression counts
22. _Scripts/combineCounts.R
23. _Scripts/consensusCluster_main.R
24. _Scripts/clusterControl_main.R
25. _Scripts/clusterJoin_main.R
22. _Scripts/Vis_main.R
26. _Scripts/DESeq2_main.R
27. _Scripts/GSEA_main.R
28. _Scripts/consensusExpression_main.R
29. _Scripts/modeltINIT_main.m																[25]

Network analysis
30. _Scripts/combineCounts_ADPD-only.R
31. _Scripts/networkAnalysis.html

Zebrafish analysis
32. _Scripts/kallisto_index_dre.sh
33. _Scripts/zeb_kallisto.sh
34. _Scripts/zeb_cts.R
35. _Scripts/zeb_DESeq2.R
36. _Scripts/zeb_GSEA.R
37. _Scripts/zeb_DESeq2_sep.R
38. _Scripts/zeb_consensusExpression.R
39. _Scripts/zeb_ENSDART2entrez.R
40. _Scripts/zeb_modeltINIT.m																[25]

================================================================================== CITATIONS ====

[ 1] McInnes L et al, 2018, "UMAP: Uniform manifold approximation  and projection." J Open Source
     Softw, 3: 381. DOI: 10.21105/joss.00861.
[ 2] Väremo L et al, 2013, “Enriching the gene set analysis  of genome-wide data by incorporating
     directionality of gene expression and combining statistical hypotheses and methods.” Nucleic
	 Acids Res, 41(8): 4378-4391. DOI: 10.1093/nar/gkt111.
[ 3] Love MI et al,  2014, “Moderated  estimation of fold  change and dispersion for RNA-seq data
     with DESeq2.” Genome Biol, 15: 550. DOI: 10.1186/s13059-014-0550-8.
[ 4] Wilkerson  DM  and Hayes ND  (2010).  “ConsensusClusterPlus:  a class  discovery  tool  with
     confidence   assessments   and   item   tracking.”   Bioinformatics,    26(12):   1572-1573.
	 DOI: 10.1093/bioinformatics/btq170
[ 5] van der Maaten LJP and Hinton GE, 2008,  "Visualizing high-dimensional data using  t-SNE", J
     Mach Learn Res, 9: 2579-2605.
[ 6] Gift  N et al,  2010,  "MetabolAnalyze:  Probabilistic  principal  components  analysis  for
	 metabolomic data." R package version 1.3.1.
[ 7] van Dijk D et al,  2018, "Recovering  gene  interactions  from  single-cell  data using data
	 diffusion," Cell, 174(3): 716-729e27. DOI: 10.1016/j.cell.2018.05.061.
[ 8] Robinson MD et al, 2010, "edgeR: a Bioconductor package for differential expression analysis
     of    digital     gene     expression     data."     Bioinformatics,     26(1):     139-140.
	 DOI: 10.1093/bioinformatics/btp616.
[ 9] Richie ME et al, 2015, "limma powers differential expression analyses for RNA-sequencing and
 	 microarray studies.” Nucleic Acids Res, 43(7): e47. DOI: 10.1093/nar/gkv007.
[10] Durinck S et al,  2009, "Mapping  identifiers for  the integration of  genomic datasets with
	 the R/Bioconductor package biomaRt." Nat Protoc, 4: 1184–1191. DOI: 10.1038/nprot.2009.97.
[11] Durinck S  et al,  2005,  "BioMart  and  Bioconductor:  a powerful  link between  biological
     databases   and    microarray    data   analysis."    Bioinformatics,   21(16):   3439-3440.
	 DOI: 10.1093/bioinformatics/bti525.
[12] Bray NL et al, 2016, "Near-optimal probablistic RNA-seq quantification." Nat Biotechnol, 34:
	 525-527. DOI: 10.1038/nbt.3519.
[13] Heirendt L et al, 2019,  "Creation and analysis of biochemical  constraint-based models: the
     COBRA Toolbox v3.0." Nat Protoc, 14: 639-702. DOI: 10.1038/s41596-018-0098-2.
[14] Wang H et al, 2018, "RAVEN 2.0: A versatile toolbox for metabolic network reconstruction and
	 a  case   study  on   Streptomyces   coelicolor."   PLoS   Comput  Biol   14(10):  e1006541.
	 DOI: 10.1371/journal.pcbi.1006541.
[15] Agren et al, 2012,  "Reconstruction of  genome-scale active  metabolic networks for 69 human
     cell   types  and  16  cancer  types   using  INIT."   PLoS  Comput  Biol,  8(5):  e1002518.
	 DOI: 10.1371/journal.pcbi.1002518.
[16] Bornstein BJ  et al,  2008,  "LibSBML:  An API  library  for  SBML."  Bioinformatics, 24(6):
	 880-881. DOI: 10.1093/bioinformatics/btn051.
[17] Marbach  D et al,  2016,  "Tissue-specific   regulatory  circuits  reveal  variable  modular
	 perturbations across complex diseases." Nat Meth, 13: 366-370, DOI: 10.1038/nmeth.3799.
[18] Luck K et al, 2020, "A reference map of the  human binary protein interactome." Nature, 580:
	 402-408. DOI: 10.1038/s41586-020-2188-x.
[19] Uhlén M et al, 2015, "Tissue-based map of the human proteome."  Science, 347(6220): 1260409.
	 DOI: 10.1126/science.1260419.
[20] Rajkumar AP et al, 2019,  "Postmortem cortical transcriptomics  of Lewy body dementia reveal
	 mitochondrial dysfunction  and lack of neuroinflammation."  Am J Geriatr Psychiatry,  28(1):
	 75-86. DOI: 10.1016/j.jagp.2019.06.007.
[20] Bennett DA  et al,  2018,  "Religious Orders  Study and  Rush Memory  and Aging Project."  J
	 Alzheimers Dis, 64(Suppl 1): S161-S189. DOI: 10.3233/JAD-179939.
[21] Zhang Y et al,  2005,  "Transcriptional  analysis  of multiple brain  regions in Parkinson's
	 disease supports  the involvement  of specific  protein  processing,  energy metabolism, and
	 signaling pathways, and suggests novel disease mechanisms."  Am J Med Genet B Neuropsychiatr
	 Genet, 137B(1): 5-16. DOI: 10.1002/ajmg.b.30195
[22] Zheng B et al,  2010,  "PGC-1α,  a potential  therapeutic  target for  early intervention in
	 Parkinson's disease." Sci Transl Med, 2(52): 52ra73. DOI: 10.1126/scitranslmed.3001059.
[23] Yates  AD   et  al,   2019,   "Ensembl   2020."   Nucleic  Acids  Res,   48(D1):  D682-D688.
     DOI: 10.1093/nar/gkz966.
[24] van Steijn et al, 2019,  "Predicting metabolism from  gene expression  in an improved whole-
	 genome   metabolic    network   model  of  Danio   rerio."    Zebrafish,   16(4):   348-362.
	 DOI: 10.1089/zeb.2018.1712.
[25] Baloni P et al, 2020,  "Metabolic network analysis  reveals altered bile  acid synthesis and
	 metabolism     in    Alzheimer's     disease."     Cell    Rep    Med,     1(8),     100138.
	 DOI: 10.1016/j.xcem.2020.100138

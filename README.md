# myodiff

Analysis and making of figures for C2C12 myodiff RNA-seq paper

<img align="center" src="https://i.imgur.com/RqxvSJf.png">

## Contacts

[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), 2016-2018, predeus at gmail dot com 

## Contents 

The repository contains processed data from the [GSE108503](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108503) RNA-Seq dataset (14 samples, 7 conditions).

It also contains all the code used in the downstream analysis (differential expression, pathway analysis). 

## Functions

* deseq_de - runs DESeq2 using expression matrix and a pair of conditions
* mass_deseq - runs many comparisons (via deseq_de) using a provided contrasts file
* pca12 - makes a well-annotated PCA plot 
* mass_fgsea_deseq - runs pre-ranked fGSEA on all outputs of DESeq2
* eliminatePathways2 - attempts to group together similar pathways among all significant results

Source the functions in *mydiff_functions.R* before running *myodiff_paper_analysis.R*.


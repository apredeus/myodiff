# Supplementary information for C2C12 myodifferentiation RNA-seq study

Analysis and scripts used in making of the manuscript figures

<img align="center" src="https://i.imgur.com/RqxvSJf.png">

## Contacts

[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), 2016-2018, predeus at gmail dot com 

## Publications

This work has been reported at Frontiers of Cardiovascular Biology conference (Vienna, April 20-22nd, 2018): 

*A. Predeus; OA. Ivanova; NV. Khromova; AM. Kiselev; DE. Polev; NA. Smolina; AA. Kostareva; RI. Dmitrieva.* [Pathway analysis of RNA-sequencing of various stages of myodifferentiaion identifies conditions favoring type I and type II fibers, and highlights increased efficiency of combined differentiation.](https://academic.oup.com/cardiovascres/article/114/suppl_1/S16/4981026)

Full manuscript is currently in preparation. 

## Funding 

This work has been funded by the Russian Science Foundation, grant # 16-15-10178.

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


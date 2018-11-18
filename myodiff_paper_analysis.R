stringsAsFactors = F
setwd("~/myodiff/data")

library(rapportools)
library(data.table)
library(DESeq2)
library(limma)
library(fgsea)
library(BiocParallel)
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(reshape)

source("~/myodiff/myodiff_functions.R")

table          <- read.table("Project_RI.sum_expression.tsv", header=T)
exp            <- table[,!names(table) %in% c("Gene_id", "Symbol", "Gene_type")]
ann            <- table[, c(1:3)]
cond           <- read.table("Conditions.txt", header=T, row.names=1)
head(table)
colSums(exp)

## PCA

dds_exp        <- DESeqDataSetFromMatrix(countData=round(exp,0),colData=cond,design = ~ Condition)
rlog_exp       <- rlog(dds_exp,blind=T)
rlog_exp       <- as.data.frame(assay(rlog_exp))  

pca12(rlog_exp,cond,ann,"Condition","Row.names")

## differential expression analysis with DEseq2

rownames(ann)  <- ann$Gene_id
ann$Gene_id    <- NULL
rownames(exp)  <- rownames(ann)

DE             <- mass_deseq(exp,cond,ann,"Contrasts.txt","Myo_de")
## this list now containts outputs of all 15 comparisons 
## written tables have all differential genes at FDR 10% cutoff

## now run pre-ranked GSEA analysis using fGSEA http://bioconductor.org/packages/release/bioc/html/fgsea.html

h.all          <- gmtPathways("~/myodiff/data/h.all.v6.0.symbols.gmt")
h_all_gsea     <- mass_fgsea_deseq(DE,h.all,1000000,"Contrasts.txt","RI_fgsea_h")

## some nice GSEA curves

rnk1              <- aggregate(stat ~ Symbol, data=DE[[1]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk1              <- setNames(rnk1$stat,toupper(rnk1$Symbol))
rnk2              <- aggregate(stat ~ Symbol, data=DE[[2]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk2              <- setNames(rnk2$stat,toupper(rnk2$Symbol))
rnk3              <- aggregate(stat ~ Symbol, data=DE[[3]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk3              <- setNames(rnk3$stat,toupper(rnk3$Symbol))
rnk4              <- aggregate(stat ~ Symbol, data=DE[[4]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk4              <- setNames(rnk4$stat,toupper(rnk4$Symbol))
rnk5              <- aggregate(stat ~ Symbol, data=DE[[5]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk5              <- setNames(rnk5$stat,toupper(rnk5$Symbol))
rnk6              <- aggregate(stat ~ Symbol, data=DE[[6]], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk6              <- setNames(rnk6$stat,toupper(rnk6$Symbol))

plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk1) + labs(title="H: Adipogenesis")
plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk2) + labs(title="H: Adipogenesis")
plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk3) + labs(title="H: Adipogenesis")
plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk4) + labs(title="H: Adipogenesis")
plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk5) + labs(title="H: Adipogenesis")
plotEnrichment(h.all[["HALLMARK_ADIPOGENESIS"]], rnk6) + labs(title="H: Adipogenesis")

plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk1) + labs(title="H: Myogenesis")
plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk2) + labs(title="H: Myogenesis")
plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk3) + labs(title="H: Myogenesis")
plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk4) + labs(title="H: Myogenesis")
plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk5) + labs(title="H: Myogenesis")
plotEnrichment(h.all[["HALLMARK_MYOGENESIS"]], rnk6) + labs(title="H: Myogenesis")

plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk1) + labs(title="H: OxPhos")
plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk2) + labs(title="H: OxPhos")
plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk3) + labs(title="H: OxPhos")
plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk4) + labs(title="H: OxPhos")
plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk5) + labs(title="H: OxPhos")
plotEnrichment(h.all[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]], rnk6) + labs(title="H: OxPhos")

plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk1) + labs(title="H: Glycolysis")
plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk2) + labs(title="H: Glycolysis")
plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk3) + labs(title="H: Glycolysis")
plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk4) + labs(title="H: Glycolysis")
plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk5) + labs(title="H: Glycolysis")
plotEnrichment(h.all[["HALLMARK_GLYCOLYSIS"]], rnk6) + labs(title="H: Glycolysis")

## diff expressed gene #s 

plot_de <- list()

plot_de[[1]]     <- read.table("Myo_de.day0_Control_vs_day2_HS.padj0.1.deseq2.tsv",header=T,row.names=1)
plot_de[[2]]     <- read.table("Myo_de.day0_Control_vs_day2_AD.padj0.1.deseq2.tsv",header=T,row.names=1)
plot_de[[3]]     <- read.table("Myo_de.day0_Control_vs_day2_HS_AD.padj0.1.deseq2.tsv",header=T,row.names=1)
plot_de[[4]]     <- read.table("Myo_de.day0_Control_vs_day7_HS.padj0.1.deseq2.tsv",header=T,row.names=1)
plot_de[[5]]     <- read.table("Myo_de.day0_Control_vs_day7_AD.padj0.1.deseq2.tsv",header=T,row.names=1)
plot_de[[6]]     <- read.table("Myo_de.day0_Control_vs_day7_HS_AD.padj0.1.deseq2.tsv",header=T,row.names=1)

name_list <- c("day2_DM1","day2_DM2","day2_DM3","day7_DM1","day7_DM2","day7_DM3")
df        <- data.frame(exp=character(),
                 type=character(), 
                 count=numeric(), 
                 stringsAsFactors=FALSE) 

for (i in 1:6) { 
  ngt_2 <- dim(plot_de[[i]][plot_de[[i]]$log2FoldChange > 2,])[1]
  ngt_4 <- dim(plot_de[[i]][plot_de[[i]]$log2FoldChange > 4,])[1]
  nlt_2 <- dim(plot_de[[i]][plot_de[[i]]$log2FoldChange < -2,])[1] ## watch that space yo
  nlt_4 <- dim(plot_de[[i]][plot_de[[i]]$log2FoldChange < -4,])[1]
  ## exclude the stronger changing genes (for stacked plots to make sense)
  ngt_2 <- ngt_2 - ngt_4
  nlt_2 <- nlt_2 - nlt_4 
  df[nrow(df) + 1,] = list(name_list[[i]],"ngt_2",ngt_2)
  df[nrow(df) + 1,] = list(name_list[[i]],"ngt_4",ngt_4)
  df[nrow(df) + 1,] = list(name_list[[i]],"nlt_2",-nlt_2)
  df[nrow(df) + 1,] = list(name_list[[i]],"nlt_4",-nlt_4)
  cat(sprintf("%s\t%d\t%d\t%d\t%d\n",name_list[[i]],ngt_2,ngt_4,nlt_2,nlt_4))
}

ggplot() + geom_bar(data = df, aes(x=exp, y=count, fill=type),stat = "identity") +
  scale_fill_brewer(type = "seq", palette = 1) + coord_flip() + scale_x_discrete(limits = rev(levels(as.factor(df$exp))))


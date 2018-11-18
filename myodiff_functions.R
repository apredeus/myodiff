pca12 = function (exp_table,cond_table,ann_table,color,label) {
  pca         <- prcomp(t(exp_table))
  print(summary(pca))
  var_PC1     <- summary(pca)$importance[2,1]
  var_PC2     <- summary(pca)$importance[2,2]
  sumvar      <- var_PC1+var_PC2
  
  loads       <- pca$rotation
  loads_ann   <- merge(loads,ann_table,by="row.names",all.x=T)
  loads_PC1   <- loads_ann[,colnames(loads_ann) %in% c("PC1","Symbol","Row.names")]
  loads_PC2   <- loads_ann[,colnames(loads_ann) %in% c("PC2","Symbol","Row.names")]
  
  loads_PC1$PC1 <- abs(loads_PC1$PC1)
  loads_PC2$PC2 <- abs(loads_PC2$PC2)
  
  loads_PC1   <- loads_PC1[order(loads_PC1$PC1,decreasing=T),]
  loads_PC2   <- loads_PC2[order(loads_PC2$PC2,decreasing=T),]
  
  #loads_PC1   <- loads_PC1[complete.cases(loads_PC1),]
  #loads_PC2   <- loads_PC2[complete.cases(loads_PC2),]
  
  top20_PC1   <- head(loads_PC1,n=20)
  top20_PC2   <- head(loads_PC2,n=20)
  
  print.data.frame(top20_PC1,row.names=F)
  print.data.frame(top20_PC2,row.names=F)
  
  pcat        <- pca$x[,1:2]
  pcat_ann    <- merge(pcat,cond_table,by="row.names")
  
  title <- paste("PC1-PC2 plot, ",sumvar*100,"% variance explained",sep="")
  xlab  <- paste("PC1, ",var_PC1*100,"% variance",sep="")
  ylab  <- paste("PC2, ",var_PC2*100,"% variance",sep="")
  
  
  ggplot(pcat_ann, aes(PC1, PC2)) + geom_point(aes(color = pcat_ann[,colnames(pcat_ann) %in% color]),size=3) + 
    geom_text_repel(data=pcat_ann,aes(label=pcat_ann[,colnames(pcat_ann) %in% label])) +
    ggtitle(title) + labs(x=xlab,y=ylab) + scale_colour_discrete(name=color)
}

mass_deseq = function (exp_table,cond_table,ann_table,pairs_file,tag) {
  library(DESeq2)
  library(data.table)
  DE       <- list()
  pairs    <- fread(pairs_file) ## data.table object with needed column, cond1, and cond2
  N        <- nrow(pairs)
  cat(sprintf("MASS-DESEQ: Processing %d following contrasts:\n\n",N))
  print(pairs)
  cat(sprintf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n"))
  
  for (i in 1:N) {
    tryCatch({
      column   <- pairs$cond_col[i]
      cond1    <- pairs$cond1[i]
      cond2    <- pairs$cond2[i]
      cat(sprintf("MASS-DESEQ: Index: %d, cond1: %s, cond2 %s\n",i,cond1,cond2))
      DE[[i]]  <- deseq_de(exp_table,cond_table,ann_table,column,cond1,cond2)
    }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
    deg        <- DE[[i]][complete.cases(DE[[i]]) & DE[[i]]$padj <= 0.1,]
    filename   <- paste(tag,".",cond1,"_vs_",cond2,".padj0.1.deseq2.tsv",sep="")
    write.table(deg,filename,sep = "\t",quote=F,row.names=F)
  }
  return(DE)
}

deseq_de = function (exp_table,cond_table,ann_table,column,cond1,cond2) {
  library(DESeq2)
  cat(sprintf("DESEQ-DE: Processing the following conditions: column is %s, cond1 is %s, cond2 is %s.\n",column,cond1,cond2))
  cond_f1    <- cond_table[,grep(column,colnames(cond_table)),drop=F]
  cond_f2    <- cond_f1[cond_f1[,1]==cond1 | cond_f1[,1]==cond2,,drop=F]
  cat(sprintf("Truncated condition table is listed below:\n"))
  cat(sprintf("---------------------------------------\n"))
  print(cond_f2)
  cat(sprintf("---------------------------------------\n"))
  exp_f      <- exp_table[,colnames(exp_table) %in% rownames(cond_f2)]
  exp_f      <- exp_f[,row.names(cond_f2)] ## reorder appropriately
  colnames(cond_f2) <- "Cond"
  dds        <- DESeqDataSetFromMatrix(countData=round(exp_f,0),colData=cond_f2,design = ~ Cond)
  deseq      <- DESeq(dds)
  res        <- results(deseq, contrast=c("Cond",cond2,cond1))
  resord     <- as.data.frame(res[order(res$padj),])
  resord     <- merge(resord, ann_table, by = "row.names", all.x = T)
  resord     <- resord[, c("Row.names","Symbol","Gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" )]
  resord     <- resord[order(resord$padj), ]
  tryCatch({
    resord$Row.names  <- colsplit(as.character(resord$Row.names),split="\\.",names=c('a','b'))$a
  }, error=function(e){cat("WARNING: No splitting done, your Ensembl IDs are already formatted correctly.\n")})
  colnames(resord)[1]   <- "Ensembl_ID"
  
  cat(sprintf("DESEQ-DE: following is the summary of differential expression:\n"))
  cat(sprintf("---------------------------------------\n"))
  summary(res)
  cat(sprintf("======================================================\n\n\n"))
  
  return(resord) 
}

mass_fgsea_deseq = function (DE,pathways,nperm,pairs_file,tag) {
  library(DESeq2)
  library(fgsea)
  library(rapportools)
  library(data.table)
  ## you can just use the letter here
  cat(sprintf("fGSEA-DESEQ: finding enrichments using %d label shuffles.\n",nperm))
  
  ## input is a list of limma DE outputs - see above 
  ## MAKE SURE YOU USE THE SAME PAIRS FILE AS IN MASS LIMMA/DESEQ2 
  GSEA     <- list()
  pairs    <- fread(pairs_file) ## data.table object with needed column, cond1, and cond2
  N        <- nrow(pairs)
  
  for (i in 1:N) {
    cond1  <- pairs$cond1[i]
    cond2  <- pairs$cond2[i]
    cat(sprintf("fGSEA with collapsing: processing tag %s, cond1 is %s, cond2 is %s\n",i,cond1,cond2))
    
    diffexp          <- DE[[i]]
    diffexp          <- diffexp[complete.cases(diffexp),]
    diffexp          <- diffexp[diffexp$Gene_type=="protein_coding",]
    
    rnk              <- aggregate(stat ~ Symbol, data=diffexp, function(x) ifelse(mean(x)>0,max(x),min(x)))
    rnk              <- setNames(rnk$stat,toupper(rnk$Symbol))
    gsea             <- fgsea(pathways, rnk, minSize=15, maxSize=500, nperm=nperm, nproc=4)
    sigpath          <- gsea[padj<0.1][order(-NES),]
    sigpath$leadingEdge    <- vapply(sigpath$leadingEdge, paste, collapse = ", ", character(1L))
    
    elim             <- eliminatePathways2(gsea,rnk,pathways)
    sigpath$Main     <- ifelse(sigpath$pathway %in% elim$pathway, "YES","-")
    setcolorder(sigpath, c(1,9,2,3,4,5,6,7,8))
    GSEA[[i]]        <- sigpath
    filename         <- paste(tag,".",cond1,"_vs_",cond2,".sigpath.gsea.tsv",sep="")
    write.table(sigpath,filename,sep = "\t",quote=F,row.names=F)
  }
  return(GSEA)
}

eliminatePathways2 = function (fgseaRes,ranks,pathways,padj.threshold=0.1) { 
  library(fgsea)
  pathways <- lapply(pathways, intersect, y=names(ranks))
  pathways <- pathways[sapply(pathways, length) >= 15]
  pathways <- pathways[sapply(pathways, length) <= 500]
  
  universe <- names(ranks)
  fgseaRes[, leadingEdge := NULL]
  
  nonInterestingReason <- setNames(rep(NA, length(pathways)), names(pathways))
  nonInterestingReason[fgseaRes[padj > padj.threshold, pathway]] <- ""
  
  i <- 2
  for (i in seq_along(pathways)) {
    p <- names(pathways)[i]
    if (!is.na(nonInterestingReason[p])) {
      next
    }
    message(p)
    
    pathwaysToCheck <- setdiff(names(which(is.na(nonInterestingReason))), p)
    messagef("checks left: <= %s", length(setdiff(pathwaysToCheck, names(head(pathways, i)))))
    
    if (length(pathwaysToCheck) == 0) {
      break
    }
    
    minPval   <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)
    u1        <- setdiff(universe, pathways[[p]])
    fgseaRes1 <- fgsea(pathways = pathways[pathwaysToCheck], stats=ranks[u1], nperm=1000, maxSize=length(u1)-1, minSize=1, nproc=4)
    minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway], fgseaRes1$pval)
    u2        <- pathways[[p]]
    fgseaRes2 <- fgsea(pathways = pathways[pathwaysToCheck], stats=ranks[u2], nperm=1000, maxSize=length(u2)-1, minSize=1, nproc=4)
    minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway], fgseaRes2$pval)
    nonInterestingReason[names(which(minPval > 0.1))] <- p
  }
  
  mainPathways <- names(which(is.na(nonInterestingReason)))
  mainTable    <- fgseaRes[match(mainPathways, pathway),]
  mainTable    <- mainTable[order(mainTable$NES),]
  return(mainTable)
}
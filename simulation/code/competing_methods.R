## Functions for simulation
## Each function returns the indices of rejections.

## In our simulation, we focus on  
## numerical or two-level categorical variable of interest. Please put the variable 
## of interest at first variable position, i.e., if the variable of interest is "u" 
## and confounder is "x", use **formula <- 'u+x'** instead of **formula <- 'x+u'**.

## We disable zero treatments of all competing methods.

##############################################
library(LinDA)
linda.fun <- function(otu.tab, meta, formula, alpha) {
  res <- linda(otu.tab, meta, paste('~', formula))
  pval <- res$output[[1]]$pvalue
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

library(ANCOMBC)
library(phyloseq)
## If use **zero_cut = 1**, the taxa with 100% zeros will still be removed.
## So use **zero_cut = 1.1** to ensure to keep all the taxa.
## Use **lib_cut = 1** to keep all the samples. 
ancombc.fun <- function(otu.tab, meta, formula, alpha) {
  OTU = otu_table(otu.tab, taxa_are_rows = TRUE)
  META = sample_data(meta)
  PHYSEQ = phyloseq(OTU, META)
  
  res <- ancombc(phyloseq = PHYSEQ, formula = formula, p_adj_method = "BH", 
                 zero_cut = 1.1, lib_cut = 1, conserve = TRUE)
  
  pval <- res$res$p_val[, 1]
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

library(ANCOMBC)
library(phyloseq)
## Under setting S6C0, we investigate the effect of zero treatment of ANCOM-BC
## Need to modify the indices.
## Error if **case == 'C1'**
## Error if some samples are removed if meta has only one column. Add a column.
ancombc.fun2 <- function(otu.tab, meta, formula, alpha) {
  OTU = otu_table(otu.tab, taxa_are_rows = TRUE)
  fake <- rep(1, nrow(meta))
  META = sample_data(cbind(meta, fake))
  PHYSEQ = phyloseq(OTU, META)
  
  res <- ancombc(phyloseq = PHYSEQ, formula = formula, p_adj_method = "BH", 
                 zero_cut = 0.9, lib_cut = 1000, struc_zero = TRUE,
                 group = 'u', conserve = TRUE)
  
  pval <- res$res$p_val[, 1]
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  if(length(pval) != nrow(otu.tab)) {
    rej <- as.numeric(gsub('taxon', '', rownames(pval)))[rej]
  }
  return(rej)
}

library(ALDEx2)
## Taxa with 100% zeros will be removed. Need to modify the indices.
aldex2.fun <- function(otu.tab, meta, formula, alpha) {
  tt <- FALSE
  if(length(all.vars(as.formula(paste('~', formula)))) == 1) {
    if(!is.numeric(meta[, formula])) {
      tt <- TRUE
    }
  }
  if(tt) {
    res <- aldex(otu.tab, meta[, formula], test = "t", effect = FALSE, denom = "all")
    pval <- res$wi.ep
  } else {
    design <- model.matrix(as.formula(paste('~', formula)), meta)
    x <- aldex.clr(otu.tab, design, denom = "all")
    res <- aldex.glm(x, design)
    pval <- res[, 8]
  }
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  if(length(pval) != nrow(otu.tab)) {
    rej <- as.numeric(gsub('taxon', '', rownames(res)))[rej]
  }
  return(rej)
}

library(DESeq2)
library(GMPR)
deseq2.fun <- function(otu.tab, meta, formula, alpha) {
  dds <- DESeqDataSetFromMatrix(countData = otu.tab,
                                colData = meta,
                                design = as.formula(paste('~', formula)))
  gmpr.size.factor <- GMPR(t(otu.tab))
  gmpr.size.factor[which(is.na(gmpr.size.factor))] <- 1
  sizeFactors(dds) <- gmpr.size.factor
  dds <- DESeq(dds)
  
  res <- results(dds, name = resultsNames(dds)[2])
  pval <- res$pvalue
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

library(edgeR)
library(GMPR)
edgeR.fun <- function(otu.tab, meta, formula, alpha) {
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  gmpr.size.factor <- GMPR(t(otu.tab))
  gmpr.size.factor[which(is.na(gmpr.size.factor))] <- 1
  data <- DGEList(counts = otu.tab, samples = meta, norm.factors = gmpr.size.factor)
  
  data <- estimateDisp(data, design)
  fit <- glmQLFit(data, design)
  qlf <- glmQLFTest(fit, coef = 2)
  pval <- qlf$table$PValue
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

library(metagenomeSeq) 
## Possible error: **Error in cumNormStatFast(data) : 
## Warning sample with one or zero features**
## metagenomeSeq currently can only analyze two-level categorical variable of interest
## with no confounders.
metagenomeSeq.fun <- function(otu.tab, meta, formula, alpha) {
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  data <- newMRexperiment(otu.tab, phenoData = AnnotatedDataFrame(meta))
  data <- cumNorm(data, p = cumNormStatFast(data))
  res <- fitFeatureModel(data, design, coef = 2)
  pval <- as.vector(res@pvalues)
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
} 

library(metagenomeSeq) 
## Possible error: **Error in cumNormStatFast(data) : 
## Warning sample with one or zero features**
## Possible error: **Error in if (max(df.residual) == 0) 
## stop("No residual degrees of freedom in linear model fits") : 
## missing value where TRUE/FALSE needed** (when there are taxa with 100% zeros)
## metagenomeSeq2 (old version) usually has severe FDR inflation
metagenomeSeq.fun2 <- function(otu.tab, meta, formula, alpha){
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  data <- newMRexperiment(otu.tab, phenoData = AnnotatedDataFrame(meta))
  data <- cumNorm(data, p = cumNormStatFast(data))
  res <- fitZig(data, design)
  pval <- as.vector(res@eb$p.value[, 2])
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
} 

library(GMPR)
## Ignore confounders (if there are) in formula,  only consider the variable of
## interest. If the variable of interest is two-level categorical, use wilcoxon;
## if the variable of interest is numerical, use spearman.
wilco.spear.fun <- function(otu.tab, meta, formula, alpha) {
  gmpr.size.factor <- GMPR(t(otu.tab))
  gmpr.size.factor[which(is.na(gmpr.size.factor))] <- 1
  Y <- t(t(otu.tab) / gmpr.size.factor)
  m <- nrow(Y)
  
  voi <- meta[, all.vars(as.formula(paste('~', formula)))[1]]
  if(!is.numeric(voi)) {
    tmp <- unique(voi)
    Y1 <- as.matrix(Y[, voi == tmp[1]])
    Y2 <- as.matrix(Y[, voi == tmp[2]])
    pval <- sapply(c(1 : m), function(i){
      wilcox.test(Y1[i, ], Y2[i, ], alternative = "two.sided")$p.value
    })
  } else {
    pval <- sapply(c(1 : m), function(i){
      cor.test(voi, Y[i, ], alternative = 'two.sided', method = 'spearman')$p.value
    })
  }
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

## Use relative abundance directly
## Ignore confounders (if there are) in formula,  only consider the variable of
## interest. If the variable of interest is two-level categorical, use wilcoxon;
## if the variable of interest is numerical, use spearman.
wilco.spear.fun2 <- function(otu.tab, meta, formula, alpha) {
  Y <- t(t(otu.tab) / colSums(otu.tab))
  m <- nrow(Y)
  
  voi <- meta[, all.vars(as.formula(paste('~', formula)))[1]]
  if(!is.numeric(voi)) {
    tmp <- unique(voi)
    Y1 <- as.matrix(Y[, voi == tmp[1]])
    Y2 <- as.matrix(Y[, voi == tmp[2]])
    pval <- sapply(c(1 : m), function(i){
      wilcox.test(Y1[i, ], Y2[i, ], alternative = "two.sided")$p.value
    })
  } else {
    pval <- sapply(c(1 : m), function(i){
      cor.test(voi, Y[i, ], alternative = 'two.sided', method = 'spearman')$p.value
    })
  }
  qval <- p.adjust(pval, method = 'BH')
  rej <- which(qval <= alpha)
  return(rej)
}

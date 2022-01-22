## Functions for real data analysis
## Each function returns a vector of p-values. 

## In our real data analysis, 
## all the datasets have binary phenotypes (CDI, IBD, RA and SMOKE).
## Please put the variable of interest at first variable position, 
## i.e., if the variable of interest is "u" 
## and confounder is "x", use **formula <- 'u+x'** instead of **formula <- 'x+u'**.

## We disable zero treatments of all competing methods.
## The input OTU table and meta data are filtered:
## 1. Taxa with prevalence (percentage of nonzeros) less than 10% are excluded. 
## 2. Samples with less than 1000 read counts are excluded.
## 3. Winsorization at quantile 0.97  

##############################################################
library(LinDA)
linda.fun <- function(otu.tab, meta, formula) {
  res <- linda(otu.tab, meta, paste('~', formula))
  pval <- res$output[[1]]$pvalue
  return(pval)
}

library(ANCOMBC)
library(phyloseq)
ancombc.fun <- function(otu.tab, meta, formula) {
  OTU = otu_table(otu.tab, taxa_are_rows = TRUE)
  META = sample_data(meta)
  PHYSEQ = phyloseq(OTU, META)
  
  res <- ancombc(phyloseq = PHYSEQ, formula = formula, p_adj_method = "BH", 
                 zero_cut = 1.1, lib_cut = 1, conserve = TRUE)
  
  pval <- res$res$p_val[, 1]
  return(pval)
}

library(ALDEx2)
aldex2.fun <- function(otu.tab, meta, formula) {
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
  return(pval)
}

library(DESeq2)
library(GMPR)
deseq2.fun <- function(otu.tab, meta, formula) {
  dds <- DESeqDataSetFromMatrix(countData = otu.tab,
                                colData = meta,
                                design = as.formula(paste('~', formula)))
  gmpr.size.factor <- GMPR(t(otu.tab))
  gmpr.size.factor[which(is.na(gmpr.size.factor))] <- 1
  sizeFactors(dds) <- gmpr.size.factor
  dds <- DESeq(dds)

  res <- results(dds, name = resultsNames(dds)[2])
  pval <- res$pvalue
  return(pval)
}

library(edgeR)
library(GMPR)
edgeR.fun <- function(otu.tab, meta, formula) {
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  gmpr.size.factor <- GMPR(t(otu.tab))
  gmpr.size.factor[which(is.na(gmpr.size.factor))] <- 1
  data <- DGEList(counts = otu.tab, samples = meta, norm.factors = gmpr.size.factor)

  data <- estimateDisp(data, design)
  fit <- glmQLFit(data, design)
  qlf <- glmQLFTest(fit, coef = 2)
  pval <- qlf$table$PValue
  return(pval)
}

library(metagenomeSeq) 
## Possible error: **Error in cumNormStatFast(data) : 
## Warning sample with one or zero features**
## metagenomeSeq currently can only analyze two-level categorical variable of interest
## with no confounders.
metagenomeSeq.fun <- function(otu.tab, meta, formula) {
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  data <- newMRexperiment(otu.tab, phenoData = AnnotatedDataFrame(meta))
  data <- cumNorm(data, p = cumNormStatFast(data))
  res <- fitFeatureModel(data, design, coef = 2)
  pval <- as.vector(res@pvalues)
  return(pval)
} 

library(metagenomeSeq)
## Possible error: **Error in cumNormStatFast(data) :
## Warning sample with one or zero features**
## Possible error: **Error in if (max(df.residual) == 0)
## stop("No residual degrees of freedom in linear model fits") :
## missing value where TRUE/FALSE needed**
## metagenomeSeq2 (old version) usually has severe FDR inflation
metagenomeSeq.fun2 <- function(otu.tab, meta, formula){
  design <- model.matrix(as.formula(paste('~', formula)), meta)
  data <- newMRexperiment(otu.tab, phenoData = AnnotatedDataFrame(meta))
  data <- cumNorm(data, p = cumNormStatFast(data))
  res <- fitZig(data, design)
  pval <- as.vector(res@eb$p.value[, 2])
  return(pval)
}
  
library(GMPR)
## Ignore confounders (if there are) in formula,  only consider the variable of
## interest. If the variable of interest is two-level categorical, use wilcoxon;
## if the variable of interest is numerical, use spearman.
wilco.spear.fun <- function(otu.tab, meta, formula) {
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
  return(pval)
}

library(Maaslin2)
maaslin2.fun <- function(otu.tab, meta, fix.eff, rand.eff) {
  res <- Maaslin2(
    input_data = otu.tab,
    input_metadata = meta,
    output = "masslin2", 
    min_prevalence = 0,
    fixed_effects = fix.eff,
    random_effects = rand.eff,
    plot_heatmap = FALSE,
    plot_scatter = FALSE)
  ind <- which(res$results$metadata == fix.eff[1])
  ord <- match(rownames(otu.tab), res$results$feature[ind])
  pval <- res$results$pval[ind][ord]
  return(pval)
}

source('competing_methods.R')

data.list <- readRDS("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/CDI_IBD_RA_SMOKE.rds")

winsor.fun <- function(Y, quan) {
  N <- colSums(Y)
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

preprocess.fun <- function(otu.tab, meta, prev.cut = 0.1, lib.cut = 1000, 
                           winsor.quan = 0.97) {
  keep.sam <- colSums(otu.tab) >= lib.cut
  Y <- otu.tab[, keep.sam]
  Z <- as.data.frame(meta[keep.sam, ])
  names(Z) <- names(meta)
  rownames(Z) <- rownames(meta)[keep.sam]
  keep.tax <- rowSums(Y > 0) / ncol(Y) >= prev.cut
  Y <- Y[keep.tax, ]
  Y <- winsor.fun(Y, winsor.quan) 
  
  return(list(Y = Y, Z = Z, keep.sam = keep.sam, keep.tax = keep.tax))
}

#######################################
formulas <- c('Disease', 'Disease+Antibiotic', 'Disease')
pval.mat.list <- list()

for(i in 1 : 3) {
  otu.tab <- data.list[[2*(i-1)+1]]
  meta <- data.list[[2*i]]
  res <- preprocess.fun(otu.tab, meta)
  Y <- res$Y
  Z <- res$Z
  
  formula <- formulas[i]
  pmat <- array(NA, c(nrow(Y), 8))
  pmat[, 1] <- linda.fun(Y, Z, formula)
  err <- try(pmat[, 2] <- ancombc.fun(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 3] <- aldex2.fun(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 4] <- deseq2.fun(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 5] <- edgeR.fun(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 6] <- metagenomeSeq.fun(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 7] <- metagenomeSeq.fun2(Y, Z, formula), silent = TRUE)
  err <- try(pmat[, 8] <- wilco.spear.fun(Y, Z, formula), silent = TRUE)
  colnames(pmat) <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'DESeq2', 'EdgeR', 
                      'MetagenomeSeq', 'MetagenomeSeq-2', 'Wilcoxon')
  rownames(pmat) <- rownames(Y)
  pval.mat.list[[i]] <- pmat
}

i <- 4
otu.tab <- data.list[[2*(i-1)+1]]
meta <- data.list[[2*i]]
res <- preprocess.fun(otu.tab, meta)
Y <- res$Y
Z <- res$Z

pmat <- array(NA, c(nrow(Y), 3))
pmat[, 1] <- linda.fun(Y, Z, formula = 'Smoke+Sex+(1|SubjectID)')
ind <- which(Z$Site == 'Left')
pmat[, 2] <- linda.fun(Y[, ind], Z[ind, ], formula = 'Smoke+Sex')
ind <- which(Z$Site == 'Right')
pmat[, 3] <- linda.fun(Y[, ind], Z[ind, ], formula = 'Smoke+Sex')
colnames(pmat) <- c('LinDA-LMM(Both)', 'LinDA-OLS(Left)', 'LinDA-OLS(Right)')
rownames(pmat) <- rownames(Y)
pval.mat.list[[i]] <- pmat
saveRDS(pval.mat.list, 'pval.mat.list.rds')
# 
# i <- 4
# otu.tab <- data.list[[2*(i-1)+1]]
# meta <- data.list[[2*i]]
# res <- preprocess.fun(otu.tab, meta)
# Y <- res$Y
# Z <- res$Z
# p1 <- linda.fun(Y, Z, formula = 'Smoke+Sex+(1|SubjectID)')
# 
# ind <- which(meta$Site == 'Left')
# otu.tab1 <- otu.tab[,ind]
# meta1 <- meta[ind,]
# res <- preprocess.fun(otu.tab1, meta1)
# Y <- res$Y
# Z <- res$Z
# p2 <- linda.fun(Y, Z, formula = 'Smoke+Sex')
# 
# ind <- which(meta$Site == 'Right')
# otu.tab2 <- otu.tab[,ind]
# meta2 <- meta[ind,]
# res <- preprocess.fun(otu.tab2, meta2)
# Y <- res$Y
# Z <- res$Z
# p3 <- linda.fun(Y, Z, formula = 'Smoke+Sex')
# 
# q1 <- p.adjust(p1, method = 'BH')
# sum(q1 <= 0.1)
# q2 <- p.adjust(p2, method = 'BH')
# sum(q2 <= 0.1)
# q3 <- p.adjust(p3, method = 'BH')
# sum(q3 <= 0.1)

pval.mat.list <- readRDS('pval.mat.list.rds')
fix.effs <- list('Disease', c('Disease','Antibiotic'), 'Disease')

for(i in 1 : 3) {
  otu.tab <- data.list[[2*(i-1)+1]]
  meta <- data.list[[2*i]]
  res <- preprocess.fun(otu.tab, meta)
  Y <- res$Y
  Z <- res$Z
  rownames(Y) <- paste0('taxon',rownames(Y))
  
  fix.eff <- fix.effs[[i]]
  pval <- maaslin2.fun(Y, Z, fix.eff, NULL)
  pmat <- pval.mat.list[[i]]
  pmat <- cbind(pmat, pval)
  colnames(pmat)[9] <- 'Maaslin2'
  pval.mat.list[[i]] <- pmat
}

i <- 4
otu.tab <- data.list[[2*(i-1)+1]]
meta <- data.list[[2*i]]
res <- preprocess.fun(otu.tab, meta)
Y <- res$Y
Z <- res$Z
rownames(Y) <- paste0('taxon',rownames(Y))

pval1 <- maaslin2.fun(Y, Z, c('Smoke','Sex'), 'SubjectID')
ind <- which(Z$Site == 'Left')
pval2 <- maaslin2.fun(Y[, ind], Z[ind, ], c('Smoke','Sex'), NULL)
ind <- which(Z$Site == 'Right')
pval3 <- maaslin2.fun(Y[, ind], Z[ind, ], c('Smoke','Sex'), NULL)
pmat <- pval.mat.list[[4]]
pmat <- cbind(pmat,pval1,pval2,pval3)
colnames(pmat)[4:6] <- c('Maaslin2(Both)', 'Maaslin2(Left)', 'Maaslin2(Right)')
pval.mat.list[[i]] <- pmat
saveRDS(pval.mat.list, 'pval.mat.list.rds')

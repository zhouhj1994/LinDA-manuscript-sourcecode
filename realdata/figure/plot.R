pval.mat.list <- readRDS("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/code_result/pval.mat.list.rds")
library(ggplot2)

num.rej.fun <- function(pval.mat, curve.fdr.cutoffs) {
  qval.mat <- sapply(1 : ncol(pval.mat), function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  isna <- colSums(is.na(qval.mat)) == nrow(pval.mat)
  
  num.rej <- t(sapply(1 : length(curve.fdr.cutoffs), function(i)
    colSums(qval.mat <= curve.fdr.cutoffs[i], na.rm = TRUE)))
  num.rej[, isna] <- NA
  return(num.rej)
}

curve.fdr.cutoffs <- seq(0.01, 0.25, 0.01)
ind <- c(1, 2, 3, 6, 9, 8)
num.rej1 <- num.rej.fun(pval.mat.list[[1]][, ind], curve.fdr.cutoffs)
num.rej2 <- num.rej.fun(pval.mat.list[[2]][, ind], curve.fdr.cutoffs)
num.rej3 <- num.rej.fun(pval.mat.list[[3]][, ind], curve.fdr.cutoffs)
n.met <- 6
cutoff <- rep(curve.fdr.cutoffs, n.met * 3)
n.cut <- length(curve.fdr.cutoffs)
method <- factor(1 : n.met, labels = c('LinDA', 'ANCOM-BC', 
                                   'ALDEx2', 'metagenomeSeq2', 'MaAsLin2', 'Wilcoxon'))
Method <- rep(rep(method, each = n.cut), 3)
dataset <- c(rep('CDI', n.cut * n.met), rep('IBD', n.cut * n.met), rep('RA', n.cut * n.met))
value <- c(as.vector(num.rej1), as.vector(num.rej2), as.vector(num.rej3))

data.plot <- cbind.data.frame(dataset, cutoff, Method, value)
plot1 <- ggplot(data.plot, aes(x = cutoff, y = value, group = Method)) +
  geom_line(aes(color = Method, linetype = Method)) +
  geom_point(aes(color = Method, shape = Method), size = 2) +
  scale_colour_manual(values = c("red","#006600","#0099CC","#FF9966", "#0000CC", "#663399"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash","longdash"))+
  scale_shape_manual(values = c(16, 17, 15, 3, 8, 7)) + 
  xlab('Target FDR Level') +
  ylab('Number of Discoveries') +
  theme_bw(base_size = 18) +
  facet_wrap(~dataset, ncol = 3, scales = 'free') +
  theme(legend.key.width = unit(1, "cm"),
        plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
pdf('num_rej_curve.pdf', width = 22, height = 8)
print(plot1)
dev.off()

################################################################## 
pval.mat.list <- readRDS("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/code_result/pval.mat.list.rds")
library(UpSetR);require(ggplot2); require(plyr); require(gridExtra); require(grid);
fun <- function(pval.mat) {
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  venn.list <- sapply(1 : n.met, function (i) 
    which(qval.mat[, i] <= cutoff))
  otu.all <- unique(names(unlist(venn.list)))
  upset.df <- data.frame(Name = otu.all)
  for(i in 1 : n.met) {
    upset.df[, i + 1] <- otu.all %in% names(venn.list[[i]]) + 0
  }
  names(upset.df)[2 : (n.met + 1)] <- method
  upset.df$Family <- as.vector(fam[match(otu.all, names(fam))])
  
  if(filename == 'CDI') {
    queries = list(list(query = elements, params = list('Family', 'Lachnospiraceae'), active = TRUE),
                   list(query = elements, params = list('Family', 'Erysipelotrichaceae'), active = TRUE))
  } else {
    queries = NULL
  }
  pdf(paste0('upset_', filename, '.pdf'), height = 8, width = 11)
  upset(upset.df, nsets = n.met, mb.ratio = c(0.5, 0.5), sets.x.label = 'Number of Discoveries',
        order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE), queries = queries,
        point.size = 3, text.scale = c(2, 2, 2, 2, 2, 2))
  dev.off()
}

cutoff <- 0.1
method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'metagenomeSeq2', 'MaAsLin2', 'Wilcoxon')
ind <- c(1, 2, 3, 6, 9, 8)
filename <- 'CDI'
pval.mat <- pval.mat.list[[1]][, ind]
load("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/CDI.RData")
fam <- data.obj$otu.name[, 'Family']


filename <- 'RA'
pval.mat <- pval.mat.list[[3]][, ind]
load("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/RA_elife.RData")
fam <- data.obj$otu.name[, 'Family']


method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'MaAsLin2', 'Wilcoxon')
ind <- c(1, 2, 3, 9, 8)
filename <- 'IBD'
pval.mat <- pval.mat.list[[2]][, ind]
load("/Users/zhouhuijuan/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/EmilyIBD.RData")
tax.df <- as.data.frame(tax)
fam <- tax.df[, 'Family']
names(fam) <- rownames(tax.df)


method <- c('LinDA-LMM(Both)', 'LinDA-OLS(Left)', 'LinDA-OLS(Right)',
            'MaAsLin2(Both)', 'MaAsLin2(Left)', 'MaAsLin2(Right)')
ind <- c(1, 2, 3, 4, 5, 6)
filename <- 'SMOKE'
pval.mat <- pval.mat.list[[4]][, ind]
load("/Users/zhouhuijuan/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/smoker_qiita_full.RData")
tax.df <- as.data.frame(smokers$tax)
fam <- tax.df[, 'Family']
names(fam) <- rownames(tax.df)

################################################################## 
pval.mat.list <- readRDS("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/code_result/pval.mat.list.rds")
library(LinDA)
library(ggplot2)
library(ggrepel)

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

formulas <- c('~Disease', '~Disease+Antibiotic', '~Disease', '~Smoke+Sex+(1|SubjectID)')
files1 <- c('plot_volcano_CDI.pdf', 'plot_volcano_IBD.pdf', 
            'plot_volcano_RA.pdf', 'plot_volcano_SMOKE.pdf')
files2 <- c('plot_lfc_CDI.pdf', 'plot_lfc_IBD.pdf', 
            'plot_lfc_RA.pdf', 'plot_lfc_SMOKE.pdf')
ind.list <- list(c(1, 2, 3, 6, 9, 8), c(1, 2, 3, 9, 8), c(1, 2, 3, 6, 9, 8), c(1, 2, 3, 4, 5, 6))
met.list <- list(c('LinDA', 'ANCOM-BC', 'ALDEx2', 'metagenomeSeq2', 'MaAsLin2', 'Wilcoxon'),
                 c('LinDA', 'ANCOM-BC', 'ALDEx2', 'MaAsLin2', 'Wilcoxon'),
                 c('LinDA', 'ANCOM-BC', 'ALDEx2', 'metagenomeSeq2', 'MaAsLin2', 'Wilcoxon'),
                 c('LinDA-LMM(Both)', 'LinDA-OLS(Left)', 'LinDA-OLS(Right)',
                   'MaAsLin2(Both)', 'MaAsLin2(Left)', 'MaAsLin2(Right)'))
titles <- c('Disease: Case v.s. DiarrhealControl', 'Disease: Crohn\'s disease v.s. Healthy',
            'Disease: HLT v.s. NORA', 'Smoke: n v.s. y')
cutoff <- 0.1
lfc.cut <- 1

plot.fun <- function(i) {
  otu.tab <- data.list[[2*(i-1)+1]]
  meta <- data.list[[2*i]]
  res <- preprocess.fun(otu.tab, meta)
  Y <- res$Y
  Z <- res$Z
  
  pval.mat <- pval.mat.list[[i]][, ind.list[[i]]]
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))
  
  linda.obj <- linda(Y, Z, formulas[i], adaptive = FALSE, imputation = FALSE)
  
  otu.tab <- linda.obj$otu.tab.use
  taxa <- rownames(otu.tab)
  bias <- linda.obj$bias
  output <- linda.obj$output
  output <- output[[1]]
  bias <- bias[1]
  lfc <- output$log2FoldChange
  lfcSE <- output$lfcSE
  padj <- output$padj
  m <- length(taxa)
  
  ## volcano plot
  leg1 <- paste0('padj>', cutoff, ' & ', 'lfc<=', lfc.cut)
  leg2 <- paste0('padj>', cutoff, ' & ', 'lfc>', lfc.cut)
  leg3 <- paste0('padj<=', cutoff, ' & ', 'lfc<=', lfc.cut)
  leg4 <- paste0('padj<=', cutoff, ' & ', 'lfc>', lfc.cut)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1 : n]
  }
  color <- gg_color_hue(3)
  
  ind1 <- padj > cutoff & abs(lfc) <= lfc.cut
  ind2 <- padj > cutoff & abs(lfc) > lfc.cut
  ind3 <- padj <= cutoff & abs(lfc) <= lfc.cut
  ind4 <- padj <= cutoff & abs(lfc) > lfc.cut
  
  leg <- rep(NA, m)
  leg[ind1] = leg1
  leg[ind2] = leg2
  leg[ind3] = leg3
  leg[ind4] = leg4
  leg <- factor(leg, levels = c(leg1, leg2, leg3, leg4))
  taxa.sig <- rep('', m)
  taxa.sig[ind3 | ind4] <- taxa[ind3 | ind4]
  
  data.volcano <- cbind.data.frame(taxa = taxa.sig, Log2FoldChange = lfc,
                                   Log10Padj = -log10(padj), leg = leg)
  plot.volcano <- ggplot(data.volcano, aes(x = Log2FoldChange, y = Log10Padj)) +
    geom_point(aes(color = leg), size = 2) +
    geom_text_repel(aes(label = taxa), max.overlaps = Inf) +
    scale_colour_manual(values = c('darkgray', color[c(2, 3, 1)])) +
    geom_hline(aes(yintercept = -log10(cutoff)), color = 'gray', linetype = 'dashed') +
    geom_vline(aes(xintercept = -lfc.cut), color = 'gray', linetype = 'dashed') +
    geom_vline(aes(xintercept = lfc.cut), color = 'gray', linetype = 'dashed') +
    ylab('-Log10Padj') +
    ggtitle(titles[i]) +
    theme_bw(base_size = 18) +
    theme(legend.title = element_blank())
  pdf(files1[i], width = 11, height = 8)
  print(plot.volcano)
  dev.off()
  
  ## effect size plot
  ind.rej.linda <- which(qval.mat[, 1] <= cutoff)
  if(i %in% c(1, 2, 3)){
    ind.rej.linda.only <- which(qval.mat[, 1] <= cutoff & 
                                  rowSums(qval.mat[,2:n.met] <= cutoff, na.rm = TRUE) == 0)
  } else if (i  == 4) {
    ind.rej.linda.only <- which(qval.mat[, 1] <= cutoff & qval.mat[, 4] > cutoff)
  }
  if(i == 1) {
    ind.rej.other.only <- which(qval.mat[, 1] > cutoff & 
                                  rowSums(qval.mat[,2:n.met] <= cutoff, na.rm = TRUE) > 0)
  } else if (i == 2 | i == 3) {
    ind.rej.other.only <- which(qval.mat[, 1] > cutoff & 
                                  rowSums(qval.mat[,2:n.met] <= cutoff, na.rm = TRUE) > 1)
  } else if(i ==  4) {
    ind.rej.other.only <- NULL
  }
  
  ind.rej <- c(ind.rej.linda, ind.rej.other.only)
  rej.color <- c(rep('black', length(ind.rej.linda)), 
                 rep('blue', length(ind.rej.other.only)))
  rej.color[which(ind.rej %in% ind.rej.linda.only)] <- 'red'
  
  n.rej <- length(ind.rej)
  taxa.rej <- taxa[ind.rej]
  taxa.rej <- factor(taxa.rej, levels = taxa.rej)
  data.plot.lfc <- cbind.data.frame(Taxa = rep(taxa.rej, 2),
                                    Log2FoldChange = c(lfc[ind.rej], lfc[ind.rej] + bias),
                                    lfcSE = c(lfcSE[ind.rej], rep(NA, n.rej)),
                                    bias = rep(c('Debiased', 'Non-debiased'), each = n.rej))
  pdf(files2[i], width = 11, height = 8)
  plot.lfc <- ggplot(data.plot.lfc, aes(x = Log2FoldChange, y = Taxa)) +
    geom_point(aes(color = bias, shape = bias), size = 3) +
    geom_errorbar(aes(xmin = Log2FoldChange - 1.96 * lfcSE,
                      xmax = Log2FoldChange + 1.96 * lfcSE), width = .2) +
    geom_vline(xintercept = 0, color = 'gray', linetype = 'dashed') +
    ggtitle(titles[i]) +
    theme_bw(base_size = 18) +
    theme(axis.text.y = element_text(colour = rej.color, size = 8),
          legend.title = element_blank())
  print(plot.lfc)
  dev.off()
}
plot.fun(1)
plot.fun(2)
plot.fun(3)
plot.fun(4)














pval.mat.list <- readRDS("pval.mat.list.rds")
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
ind <- c(1, 2, 3, 6, 8)
num.rej1 <- num.rej.fun(pval.mat.list[[1]][, ind], curve.fdr.cutoffs)
num.rej2 <- num.rej.fun(pval.mat.list[[2]][, ind], curve.fdr.cutoffs)
num.rej3 <- num.rej.fun(pval.mat.list[[3]][, ind], curve.fdr.cutoffs)

cutoff <- rep(curve.fdr.cutoffs, 5 * 3)
n.cut <- length(curve.fdr.cutoffs)
method <- factor(1 : 5, labels = c('LinDA', 'ANCOM-BC', 
                                   'ALDEx2', 'MetagenomeSeq', 'Wilcoxon'))
Method <- rep(rep(method, each = n.cut), 3)
dataset <- c(rep('CDI', n.cut * 5), rep('IBD', n.cut * 5), rep('RA', n.cut * 5))
value <- c(as.vector(num.rej1), as.vector(num.rej2), as.vector(num.rej3))

data.plot <- cbind.data.frame(dataset, cutoff, Method, value)
plot1 <- ggplot(data.plot, aes(x = cutoff, y = value, group = Method)) +
  geom_line(aes(color = Method, linetype = Method)) +
  geom_point(aes(color = Method, shape = Method), size = 2) +
  scale_colour_manual(values = c("red","#006600","#0099CC","#FF9966", "#663399"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash"))+
  scale_shape_manual(values = c(16, 17, 15, 3, 7)) + 
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
pval.mat.list <- readRDS("pval.mat.list.rds")
require(VennDiagram)

fun <- function(pval.mat, cutoff, method, filename) {
  n.met <- ncol(pval.mat)
  qval.mat <- sapply(1 : n.met, function(i) 
    p.adjust(pval.mat[, i], method = 'BH'))

  venn.list <- sapply(1 : n.met, function (i) 
    which(qval.mat[, i] <= cutoff))
  
  if(n.met == 4) {
    cat.just <- list(c(-0.2,2) , c(0.7,2.5) , c(0,2.2) , c(0.2,2.4))
  } else if(n.met == 5) {
    cat.just <- list(c(0.6,2) , c(-0.5,-2) , c(-1,-1) , c(1,-1) , c(2,-2))
  } else if(n.met == 3) {
    cat.just <- list(c(1.2,-6), c(0.2,1), c(0.4,-3))
  }
  
  venn.plot <- venn.diagram(venn.list, NULL,
                            fill = c("blue", "red", "cyan", "orange", 
                                     "green", "purple")[1 : n.met], 
                            alpha = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)[1 : n.met], cex = 3, 
                            category.names = method, 
                            cat.cex = 2,
                            cat.fontface = 1, 
                            cat.default.pos = "text",
                            cat.just = cat.just)
  pdf(paste0('venn_', filename, '.pdf'), height = 8, width = 11)
  grid.newpage()
  grid.draw(venn.plot)
  dev.off()
}

cutoff <- 0.1
method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'MetagenomeSeq', 'Wilcoxon')
ind <- c(1, 2, 3, 6, 8)
filename <- 'CDI'
pval.mat <- pval.mat.list[[1]][, ind]
fun(pval.mat, cutoff, method, filename)

filename <- 'RA'
pval.mat <- pval.mat.list[[3]][, ind]
fun(pval.mat, cutoff, method, filename)

method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'Wilcoxon')
ind <- c(1, 2, 3, 8)
filename <- 'IBD'
pval.mat <- pval.mat.list[[2]][, ind]
fun(pval.mat, cutoff, method, filename)

method <- c('LinDA-LMM(Both)', 'LinDA-OLS(Left)', 'LinDA-OLS(Right)')
ind <- c(1, 2, 3)
filename <- 'SMOKE'
pval.mat <- pval.mat.list[[4]][, ind]
fun(pval.mat, cutoff, method, filename)

##################################################################
pval.mat.list <- readRDS("pval.mat.list.rds")
pval.mat <- pval.mat.list[[1]][, c(1, 2, 3, 6, 8)]
n.met <- ncol(pval.mat)
qval.mat <- sapply(1 : n.met, function(i) 
  p.adjust(pval.mat[, i], method = 'BH'))
rej.list <- sapply(1 : n.met, function (i) 
  rownames(pval.mat)[which(qval.mat[, i] <= 0.1)])

load("/Users/zhouhuijuan/Documents/diff_abundance_analysis/DAA20210315/LinDA/realdata/data/CDI.RData")
fam <- data.obj$otu.name[, 'Family']
fam[which(rownames(data.obj$otu.name) %in% rej.list[[1]])]
fam[which(rownames(data.obj$otu.name) %in% rej.list[[2]])]
fam[which(rownames(data.obj$otu.name) %in% rej.list[[3]])]
fam[which(rownames(data.obj$otu.name) %in% rej.list[[4]])]
fam[which(rownames(data.obj$otu.name) %in% rej.list[[5]])]

##################################################################
library(LinDA)

data.list <- readRDS("/Users/zhouhuijuan/Documents/diff_abundance_analysis/DAA20210315/LinDA/realdata/data/CDI_IBD_RA_SMOKE.rds")

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
directories <- c('figures_CDI', 'figures_IBD', 'figures_RA', 'figures_SMOKE')

for(i in 1 : 4) {
  otu.tab <- data.list[[2*(i-1)+1]]
  meta <- data.list[[2*i]]
  res <- preprocess.fun(otu.tab, meta)
  Y <- res$Y
  Z <- res$Z
  
  formula <- formulas[i]
  linda.obj <- linda(Y, Z, formula)
  
  rej <- rownames(linda.obj$output[[1]])[which(linda.obj$output[[1]]$reject)]
  nrej <- rownames(linda.obj$output[[1]])[-which(linda.obj$output[[1]]$reject)]
  taxa.plot <- c(rej[1], nrej[1])
  if(i == 3) {
    height <- 11
    width <- 15
  } else {
    height <- 8
    width <- 11
  }
  plot.obj <- plot.fun(linda.obj, alpha = 0.1, lfc.cut = 1, 
                       taxa.plot = taxa.plot, legend = TRUE, directory = directories[i],
                       width = width, height = height)
}

################################################################## Add
library(UpSetR);require(ggplot2); require(plyr); require(gridExtra); require(grid);
pval.mat.list <- readRDS("pval.mat.list.rds")
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
        order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE), queries = queries)
  dev.off()
}

cutoff <- 0.1
method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'MetagenomeSeq', 'Wilcoxon')
ind <- c(1, 2, 3, 6, 8)
filename <- 'CDI'
pval.mat <- pval.mat.list[[1]][, ind]
load("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/CDI.RData")
fam <- data.obj$otu.name[, 'Family']


filename <- 'RA'
pval.mat <- pval.mat.list[[3]][, ind]
load("~/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/RA_elife.RData")
fam <- data.obj$otu.name[, 'Family']


method <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'Wilcoxon')
ind <- c(1, 2, 3, 8)
filename <- 'IBD'
pval.mat <- pval.mat.list[[2]][, ind]
load("/Users/zhouhuijuan/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/EmilyIBD.RData")
tax.df <- as.data.frame(tax)
fam <- tax.df[, 'Family']
names(fam) <- rownames(tax.df)


method <- c('LinDA-LMM(Both)', 'LinDA-OLS(Left)', 'LinDA-OLS(Right)')
ind <- c(1, 2, 3)
filename <- 'SMOKE'
pval.mat <- pval.mat.list[[4]][, ind]
load("/Users/zhouhuijuan/Documents/diff_abundance_analysis/LinDA-manuscript/LinDA-result/realdata/data/smoker_qiita_full.RData")
tax.df <- as.data.frame(smokers$tax)
fam <- tax.df[, 'Family']
names(fam) <- rownames(tax.df)


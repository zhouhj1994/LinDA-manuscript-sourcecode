fdr.fun <- function(input1, input2, compr.methods){
  value <- as.vector(t(input1))
  sd <- as.vector(t(input2))
  
  sample.size.0 <- c(1, 2)
  sig.density.0 <- c(1, 2)
  sig.strength.0 <- c(1, 2, 3, 4, 5, 6)
  
  s1 <- 2; s2 <- 2; s3 <- 6
  
  sample.size <- rep(sample.size.0, each = s2 * s3)
  sig.density <- rep(rep(sig.density.0, each = s3), s1)
  sig.strength <- rep(sig.strength.0, s1 * s2)
  
  s <- s1 * s2 * s3 
  
  setting.0 <- cbind(sample.size, sig.density, sig.strength)
  
  n.met <- length(compr.methods)
  setting <- setting.0[rep(c(1 : s), each = n.met), ] 
  Method <- rep(compr.methods, s)
  data <- cbind.data.frame(setting, Method, value, sd)
  
  plot1 <- ggplot(data, aes(x = sig.strength, y = value, group = Method)) + 
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method, shape = Method), size = 2) +
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(0)) +
    facet_grid(factor(sig.density, labels = c("Sparse Signal", "Dense Signal"))
               ~factor(sample.size, labels = c("n = 50", "n = 200"))) +
    xlab("Signal Strength") +
    ylab("False Discovery Proportion") +
    geom_hline(aes(yintercept = 0.05), color = "black", linetype = "dashed") +
    theme_bw(base_size = 18) +
    ggtitle('A') +
    theme(plot.title = element_text(hjust = -0.05, vjust = -6), 
          legend.key.width = unit(1, "cm"),
          plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
  return(plot1)
}

power.fun <- function(input1, input2, compr.methods){
  value <- as.vector(t(input1))
  sd <- as.vector(t(input2))
  
  sample.size.0 <- c(1, 2)
  sig.density.0 <- c(1, 2)
  sig.strength.0 <- c(1, 2, 3, 4, 5, 6)
  
  s1 <- 2; s2 <- 2; s3 <- 6
  
  sample.size <- rep(sample.size.0, each = s2 * s3)
  sig.density <- rep(rep(sig.density.0, each = s3), s1)
  sig.strength <- rep(sig.strength.0, s1 * s2)
  
  s <- s1 * s2 * s3 
  
  setting.0 <- cbind(sample.size, sig.density, sig.strength)
  
  n.met <- length(compr.methods)
  setting <- setting.0[rep(c(1 : s), each = n.met), ] 
  Method <- rep(compr.methods, s)
  data <- cbind.data.frame(setting, Method, value, sd)
  
  plot1 <- ggplot(data, aes(x = sig.strength, y = value, group = Method)) + 
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method, shape = Method), size = 2) +
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(0)) +
    facet_grid(factor(sig.density, labels = c("Sparse Signal", "Dense Signal"))
               ~factor(sample.size, labels = c("n = 50", "n = 200"))) +
    xlab("Signal Strength") +
    ylab("True Positive Rate") +
    theme_bw(base_size = 18) +
    ggtitle('B') +
    theme(plot.title = element_text(hjust = -0.05, vjust = -6), 
          legend.key.width = unit(1, "cm"),
          plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
  return(plot1)
}


#########################################################
library(ggplot2)

fun <- function(setup) {
  n.met <- ncol(output)/4
  fdr <- output[, 1 : n.met][, ind]
  power <- output[, (n.met+1) : (2*n.met)][, ind]
  fdr.sd <- output[, (2*n.met+1) : (3*n.met)][, ind]
  power.sd <- output[, (3*n.met+1) : (4*n.met)][, ind]
  
  fdr.sd.1 <- fdr.sd * 1.96 
  fdr.sd.1[, 1 : ncol(fdr)] <- NA
  power.sd.1 <- power.sd * 1.96 
  power.sd.1[, 1 : ncol(fdr)] <- NA
  
  # power.fun(power, power.sd.1, compr.methods)
  # fdr.fun(fdr, fdr.sd.1, compr.methods)
  
  pdf(file = paste0(setup, "_power.pdf"), width = 11, height = 8)
  print(power.fun(power, power.sd.1, compr.methods))
  dev.off()
  pdf(file = paste0(setup, "_fdr.pdf"), width = 11, height = 8)
  print(fdr.fun(fdr, fdr.sd.1, compr.methods))
  dev.off()
}

compr.methods <- c("DESeq2", "EdgeR", "MetagenomeSeq-2")
ind <- c(5, 6, 8)
output <- read.table("output_S0C0.txt")
fun('S0C0others')

compr.methods <- c("Adaptive", "Pseudo-count", "Imputation")
compr.methods <- factor(compr.methods, levels = compr.methods)
ind <- 1:3
output <- read.table("output_S0C0LinDA.txt")
fun('S0C0LinDA')
output <- read.table("output_S6C0LinDA.txt")
fun('S6C0LinDA')

compr.methods <- c("ANCOM-BC-1", "ANCOM-BC-2")
ind <- 2:3
output <- read.table("output_S6C0.txt")
fun('S6C0ANCOMBC')

compr.methods <- c("LinDA-LMM", "LinDA-OLS", "CLR+LMM", "CLR+OLS")
compr.methods <- factor(compr.methods, levels = compr.methods)
ind <- c(3, 1, 4, 2)
output <- read.table("output_S7.1C0.txt")
fun('S71C0')
output <- read.table("output_S7.2C0.txt")
fun('S72C0')



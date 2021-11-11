fdr.fun <- function(input1, input2, compr.methods, 
                    color.manual, linetype.manual, shape.manual, sample.size.label) {
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
  method.0 <- factor(1 : n.met, labels = compr.methods)
  Method <- rep(method.0, s)
  data <- cbind.data.frame(setting, Method, value, sd)
  
  plot1 <- ggplot(data, aes(x = sig.strength, y = value, group = Method)) + 
    geom_line(aes(color = Method, linetype = Method)) +
    geom_point(aes(color = Method, shape = Method), size = 2) +
    scale_colour_manual(values = color.manual)+
    scale_linetype_manual(values = linetype.manual)+
    scale_shape_manual(values = shape.manual) + 
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(0)) +
    facet_grid(factor(sig.density, labels = c("Sparse Signal", "Dense Signal"))
               ~factor(sample.size, labels = sample.size.label)) +
    xlab("Signal Strength") +
    ylab("Empirical False Discovery Rate") +
    geom_hline(aes(yintercept = 0.05), color = "black", linetype = "dashed") +
    theme_bw(base_size = 18) +
    ggtitle('A') +
    theme(plot.title = element_text(hjust = -0.05, vjust = -6), 
          legend.key.width = unit(1, "cm"),
          plot.margin=unit(c(1, 1, 1, 1.5),"cm"))
  return(plot1)
}

power.fun <- function(input1, input2, compr.methods, 
                      color.manual, linetype.manual, shape.manual, sample.size.lable) {
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
  method.0 <- factor(1 : n.met, labels = compr.methods)
  Method <- rep(method.0, s)
  data <- cbind.data.frame(setting, Method, value, sd)
  
  plot1 <- ggplot(data, aes(x = sig.strength, y = value, group = Method)) + 
    geom_line(aes(color = Method, linetype = Method)) +
    geom_point(aes(color = Method, shape = Method), size = 2) +
    scale_colour_manual(values = color.manual)+
    scale_linetype_manual(values = linetype.manual)+
    scale_shape_manual(values = shape.manual) + 
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .2,
                  position = position_dodge(0)) +
    facet_grid(factor(sig.density, labels = c("Sparse Signal", "Dense Signal"))
               ~factor(sample.size, labels = sample.size.label)) +
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
  fdr <- output[, 1 : 10][, ind]
  power <- output[, 11 : 20][, ind]
  fdr.sd <- output[, 21 : 30][, ind]
  power.sd <- output[, 31 : 40][, ind]
  
  fdr.sd.1 <- fdr.sd * 1.96 
  fdr.sd.1[, 2 : ncol(fdr)] <- NA
  power.sd.1 <- power.sd * 1.96 
  power.sd.1[, 1 : ncol(fdr)] <- NA
  
  # power.fun(power, power.sd.1, compr.methods, 
  #           color.manual, linetype.manual, shape.manual, sample.size.label)
  # fdr.fun(fdr, fdr.sd.1, compr.methods, 
  #         color.manual, linetype.manual, shape.manual, sample.size.label)
  
  pdf(file = paste0(setup, "_power.pdf"), width = 11, height = 8)
  print(power.fun(power, power.sd.1, compr.methods, 
                  color.manual, linetype.manual, shape.manual, sample.size.label))
  dev.off()
  pdf(file = paste0(setup, "_fdr.pdf"), width = 11, height = 8)
  print(fdr.fun(fdr, fdr.sd.1, compr.methods, 
                color.manual, linetype.manual, shape.manual, sample.size.label))
  dev.off()
}

fun <- function(setup) {
  fdr <- output[, 1 : 10][, ind]
  power <- output[, 11 : 20][, ind]
  fdr.sd <- output[, 21 : 30][, ind]
  power.sd <- output[, 31 : 40][, ind]
  
  fdr.sd.1 <- fdr.sd * 1.96 
  fdr.sd.1[, 2 : ncol(fdr)] <- NA
  power.sd.1 <- power.sd * 1.96 
  power.sd.1[, 1 : ncol(fdr)] <- NA
  
  # power.fun(power, power.sd.1, compr.methods, 
  #           color.manual, linetype.manual, shape.manual, sample.size.label)
  # fdr.fun(fdr, fdr.sd.1, compr.methods, 
  #         color.manual, linetype.manual, shape.manual, sample.size.label)
  
  pdf(file = paste0(setup, "_power_new.pdf"), width = 11, height = 8)
  print(power.fun(power, power.sd.1, compr.methods, 
                   color.manual, linetype.manual, shape.manual, sample.size.label))
  dev.off()
  pdf(file = paste0(setup, "_fdr_new.pdf"), width = 11, height = 8)
  print(fdr.fun(fdr, fdr.sd.1, compr.methods, 
                color.manual, linetype.manual, shape.manual, sample.size.label))
  dev.off()
}

compr.methods <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'metagenomeSeq', 'Wilcoxon')
color.manual <- c("red","#006600","#0099CC","#FF9966", "#663399")
linetype.manual <- c("solid", "dashed", "dotted", "dotdash", "longdash")
shape.manual <- c(16, 17, 15, 3, 7)
ind <- c(1, 2, 4, 7, 9)
sample.size.label <- c("n = 50", "n = 200")

output <- read.table("output_S0C0.txt")
fun('S0C0')
output <- read.table("output_S1C0.txt")
fun('S1C0')
output <- read.table("output_S2C0.txt")
fun('S2C0')
output <- read.table("output_S3C0.txt")
fun('S3C0')
output <- read.table("output_S4C0.txt")
fun('S4C0')
output <- read.table("output_S6C0.txt")
fun('S6C0')
output <- read.table("output_S0C0Strong.txt")
fun('S0C0Strong')

sample.size.label <- c("n = 20", "n = 30")
output <- read.table("output_S5C0.txt")
fun('S5C0')

compr.methods <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'Spearman')
color.manual <- c("red","#006600","#0099CC", "#663399")
linetype.manual <- c("solid", "dashed", "dotted", "longdash")
shape.manual <- c(16, 17, 15, 7)
ind <- c(1, 2, 4, 9)
sample.size.label <- c("n = 50", "n = 200")
output <- read.table("output_S0C1.txt")
fun('S0C1')

compr.methods <- c('LinDA', 'ANCOM-BC', 'ALDEx2', 'Wilcoxon')
output <- read.table("output_S0C2.txt")
fun('S0C2')



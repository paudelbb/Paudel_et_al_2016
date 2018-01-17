# ===========================================================================================================
# Supp. Figure 4B: SKMEL5 Subclones Mix Model Fitting and MCMC Algorithms
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)

setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
source("~//Paudel_et_al_2016//Codes//summarySE.R") # Summarize function to get the summary statistics

# ===========================================================================================================
#                                   Read experimental data
# ===========================================================================================================
d1 <- read.csv("MCMC_ensembles_A375.csv", header=T, sep=",")
d2 <- read.csv("MCMC_ensembles_WM793.csv", header = T, sep = ",")
d3 <- read.csv("MCMC_ensembles_SKMEL5_Parental.csv", header=T, sep=",")
d4 <- read.csv("MCMC_ensembles_SKMEL19.csv", header=T, sep=",")
d5 <- read.csv("MCMC_ensembles_SKMEL28.csv" , header=T, sep=",")
d6 <- read.csv("MCMC_ensembles_WM164.csv", header=T, sep=",")
d7 <- read.csv("MCMC_ensembles_WM88.csv", header = T, sep = ",")

d1$cell = "A375"
d2$cell = "WM793"
d3$cell = "SKMEL5"
d4$cell = "SKMEL19"
d5$cell = "SKMEL28"
d6$cell = "WM164"
d7$cell = "WM88"
# ===========================================================================================================
# ===========================================================================================================
data <- rbind(d1, d2, d3, d4, d5, d6, d7)
data$fA = data$A / (data$A + data$B + data$C)
data$fB = data$B / (data$A + data$B + data$C)
data$fC = data$C / (data$A + data$B + data$C)
data1 <- subset(data, data$time ==0) # initial fraction of cells in different cell lines
data2 <- subset(data, data$time ==200) # fraction of cells in idling state in different cell lines
# ===========================================================================================================
# Plot the boxplot of the fraction of cells in idling cells
# ===========================================================================================================

data2$cell <-factor(data2$cell, levels=c("SKMEL5", "A375", "WM793", "SKMEL19", "SKMEL28", "WM164", "WM88"))
labels = unique(data2$cell)
pdf(paste0(output, "/Proportion_of_idling_cells_in_stateA_for_all.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fA~data2$cell, data=data2, ylim=c(0,1.), main="R(idling)", notch=T, xlab="", xaxt="n")
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste0(output, "/Proportion_of_idling_cells_in_stateB_for_all.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fB~data2$cell, data=data2, ylim=c(0,1.), main="S(idling)", notch=T, xaxt="n", xlab="")
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()

pdf(paste0(output, "/Proportion_of_idling_cells_in_stateC_for_all.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fC~data2$cell, data=data2, ylim=c(0,1.), main="E(idling)", notch=T, xaxt="n", xlab="")
text(x =  seq_along(labels), y = par("usr")[3]-0.05, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
dev.off()
# ===========================================================================================================




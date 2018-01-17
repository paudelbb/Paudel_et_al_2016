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
d1 <- read.csv("MCMC_ensembles_SKMEL5 Subclone01.csv", header=T, sep=",")
d2 <- read.csv("MCMC_ensembles_SKMEL5 Subclone07.csv", header=T, sep=",")
d3 <- read.csv("MCMC_ensembles_SKMEL5 Subclone10.csv", header=T, sep=",")
d4 <- read.csv("MCMC_ensembles_SKMEL5_Parental.csv", header=T, sep=",")

d1$cell = "SC01"
d2$cell = "SC07"
d3$cell = "SC10"
d4$cell = "SKMEL5"
# ===========================================================================================================
# ===========================================================================================================
data <- rbind(d1, d2, d3, d4)
data$fA = data$A / (data$A + data$B + data$C)
data$fB = data$B / (data$A + data$B + data$C)
data$fC = data$C / (data$A + data$B + data$C)
data1 <- subset(data, data$time ==0) # initial fraction of cells in different cell lines
data2 <- subset(data, data$time ==280) # fraction of cells in idling state in different cell lines

# ===========================================================================================================
# Plot the boxplot of the fraction of cells in idling cells
# ===========================================================================================================
pdf(paste0(output, "/Proportion_of_idling_cells_in_stateA.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fA~cell, data=data2, ylim=c(0,1.), main="Proportion A", notch=T)
dev.off()

pdf(paste0(output, "/Proportion_of_idling_cells_in_stateB.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fB~cell, data=data2, ylim=c(0,1.), main="Proportion B", notch=T)
dev.off()

pdf(paste0(output, "/Proportion_of_idling_cells_in_stateC.pdf"), width=3.2, height=3.2)
par(ps = 8, cex = 1, cex.axis = 1)
boxplot(fC~cell, data=data2, ylim=c(0,1.), main="Proportion C", notch=T)
dev.off()
# ===========================================================================================================
# ===========================================================================================================
data2$ratio = data2$A/data2$C
boxplot(ratio~cell, data=data2, outline=T, notch=T)
abline(h=0.273, lty=2)






# ===========================================================================================================
# Figure 4E: Model proportions compared with states categorization cFP Assay
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)
require(Hmisc)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
source("~//Paudel_et_al_2016//Codes//summarySE.R") # Summarize function to get the summary statistics
# ===========================================================================================================
# Read experimental data
# ===========================================================================================================
d1 <- read.csv("State_proportions_cFP_assay.csv", header=T, sep=",")
data <- read.csv("MCMC_model_predictions_state_proportions.csv", header=T, sep=",")

sum1 <- summarySE(d1, measurevar = "A0", groupvars = c("cell"))
sum2 <- summarySE(d1, measurevar = "B0", groupvars = c("cell"))
sum3 <- summarySE(d1, measurevar = "C0", groupvars = c("cell"))
# ===========================================================================================================
# ===========================================================================================================

t <- data[(data$cell %in% c("SC01","SC07","SC10","SKMEL5")),]
t$cell <- factor(t$cell)
t1 <- subset(data, data$cell =="SKMEL19")
t1$cell <- "WM70"
t1$cell <- factor(t1$cell)
t2 <- subset(data, data$cell =="WM88")
t2$cell <- factor(t2$cell)
data <- rbind(t, t1, t2)

# ===========================================================================================================
# Define the labels
# ===========================================================================================================

group = c("SC01","SC07","SC10","SKMEL5", "SKMEL19","WM88")
group <- noquote(group)
labels <- group

x = c(1,2,3,4,5,6)
# ===========================================================================================================
# ===========================================================================================================
pdf(paste0(output, "/Proportions of cells in stateA.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(A0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
points(x, sum1$A0, col="red", pch=8, lwd=1.5)
dev.off()

pdf(paste0(output, "/Proportions of cells in stateB.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(B0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
points(x, sum2$B0, col="red", pch=8, lwd=1.5)
dev.off()

pdf(paste0(output, "/Proportions of cells in stateC.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(C0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
points(x, sum3$C0, col="red", pch=8, lwd=1.5)
dev.off()
# ===========================================================================================================
# Additional figures if reviewer requests (m +- 0.10 *m) on state categorization
# ===========================================================================================================
output = paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "/Figures/AdditionalFigures")

pdf(paste0(output, "/Proportions of cells in stateA.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(A0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
errbar(x, sum1$A0, sum1$A0-sum1$sd, sum1$A0+sum1$sd, col="red", pch=20, lwd=1.0, errbar.col = "red", add=T)
dev.off()

pdf(paste0(output, "/Proportions of cells in stateB.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(B0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
errbar(x, sum2$B0, sum2$B0-sum2$sd, sum2$B0+sum2$sd, col="red", pch=20, lwd=1.0, errbar.col = "red", add=T)
dev.off()

pdf(paste0(output, "/Proportions of cells in stateC.pdf"), width=3.0, height=3.2)
par(ps = 10, cex = 1, cex.axis = 1)
boxplot(C0~cell, data=data, outline=T, ylim=c(0,1), notch=T, xaxt="n", col=grey.colors(6))
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
errbar(x, sum3$C0, sum3$C0-sum3$sd, sum3$C0+sum3$sd, col="red", pch=20, lwd=1.0, errbar.col = "red", add=T)
dev.off()

# ===========================================================================================================
# ===========================================================================================================

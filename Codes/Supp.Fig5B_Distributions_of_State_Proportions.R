# ===========================================================================================================
# Supp. Figure 5B: SKMEL5 distribution of rate constants
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
fig = "Supp.Fig5/State_Proportions" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)

# ===========================================================================================================
# 
# ===========================================================================================================
load("MCMC1_A375_150000_Run.RData") # A375
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data1 <- (data[,c(5,6)])
data1$cell = "A375"

load("MCMC1_SKMEL5_Parental_150000_Run.RData") # SKMEL5
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data2 <- (data[,c(5,6)])
data2$cell = "SKMEL5"

load("MCMC1_SKMEL19_150000_Run.RData") # SKMEL19
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data3 <- (data[,c(5,6)])
data3$cell = "SKMEL19"

load("MCMC1_SKMEL28_2e+05_Run.RData") # SKMEL28
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/1.6),]
data <- data.frame(MCMC1$pars)
data4 <- (data[,c(5,6)])
data4$cell = "SKMEL28"

load("MCMC1_WM88_150000_Run.RData") # WM88
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data5 <- (data[,c(5,6)])
data5$cell = "WM88"

load("MCMC1_WM164_1e_05_Run.RData") # WM164
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data6 <- (data[,c(5,6)])
data6$cell = "WM164"

load("MCMC1_WM793_150000_Run.RData") # WM793
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data7 <- (data[,c(5,6)])
data7$cell = "WM793"

load("MCMC1_A2058_150000_Run.RData") # A2058
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data8 <- (data[,c(5,6)])
data8$cell = "A2058"

data = rbind(data1, data2, data3, data4, data5, data6, data7, data8)
data$C0 = 1 - (data$A0 + data$B0)

# ===========================================================================================================
# Rank cells based on their initial proportions
# ===========================================================================================================

source("~/Paudel_et_al_2016/Codes/summarySE.R")
dfA <- summarySE(data, measurevar="A0", groupvars=c("cell"))
dfB <- summarySE(data, measurevar="B0", groupvars=c("cell"))
dfC <- summarySE(data, measurevar="C0", groupvars=c("cell"))
# 
par(mfrow=c(2,2))
mybarcol = "gray20"
par(ps = 10, cex = 1, cex.axis = 1)
mp = barplot2(dfA$A0, beside = T, 
              col=rainbow(nrow(dfA)), ylab="", 
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim=c(0, 1), 
              plot.ci = T, ci.l = dfA$A0-dfA$ci, ci.u = dfA$A0+dfA$ci, 
              plot.grid = T)
group = c(as.character((dfA$cell)))
group = noquote(group)
labels = paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.45)

mp = barplot2(dfB$B0, beside = T, 
              col=rainbow(nrow(dfA)), ylab="", 
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim=c(0, 1), 
              plot.ci = T, ci.l = dfB$B0-dfB$ci, ci.u = dfB$B0+dfB$ci, 
              plot.grid = T)
group = c(as.character((dfB$cell)))
group = noquote(group)
labels = paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.45)

mp = barplot2(dfC$C0, beside = T, 
              col=rainbow(nrow(dfA)), ylab="", 
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim=c(0, 1), 
              plot.ci = T, ci.l = dfC$C0-dfC$ci, ci.u = dfC$C0+dfC$ci, 
              plot.grid = T)
group = c(as.character((dfC$cell)))
group = noquote(group)
labels = paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.45)



# ===========================================================================================================
# Plot the distributions states A0
# ===========================================================================================================
par(ps = 12, cex = 1, cex.axis = 1)
# A0
pdf(paste0(output, "/", "A375_A0.pdf"), width=3.5, height=3.5)
hist(data1$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL5_A0.pdf"), width=3.5, height=3.5)
hist(data2$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL19_A0.pdf"), width=3.5, height=3.5)
hist(data3$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL28_A0.pdf"), width=3.5, height=3.5)
hist(data4$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM88_A0.pdf"), width=3.5, height=3.5)
hist(data5$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM164_A0.pdf"), width=3.5, height=3.5)
hist(data6$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM793_A0.pdf"), width=3.5, height=3.5)
hist(data7$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "A2058_A0.pdf"), width=3.5, height=3.5)
hist(data8$A0, breaks=seq(0, 1, by=0.005), main="A0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()


# ===========================================================================================================
# Plot the distributions of B0
# ===========================================================================================================

# B0
pdf(paste0(output, "/", "A375_B0.pdf"), width=3.5, height=3.5)
hist(data1$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL5_B0.pdf"), width=3.5, height=3.5)
hist(data2$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL19_B0.pdf"), width=3.5, height=3.5)
hist(data3$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "SKMEL28_B0.pdf"), width=3.5, height=3.5)
hist(data4$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM88_B0.pdf"), width=3.5, height=3.5)
hist(data5$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM164_B0.pdf"), width=3.5, height=3.5)
hist(data6$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "WM793_B0.pdf"), width=3.5, height=3.5)
hist(data7$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()

pdf(paste0(output, "/", "A2058_B0.pdf"), width=3.5, height=3.5)
hist(data8$B0, breaks=seq(0, 1, by=0.005), main="B0", freq = F, xlab="", xlim=c(0,1), ylim=c(0,12))
dev.off()
# ===========================================================================================================
# ===========================================================================================================

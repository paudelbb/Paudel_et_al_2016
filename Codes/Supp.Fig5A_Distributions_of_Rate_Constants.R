# ===========================================================================================================
# Supp. Figure 5A: SKMEL5 distribution of rate constants
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
fig = "Supp.Fig5/Rate_Constants" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)

files1 = grep("prior", dir(), value = T)
a2058 = read.csv("priorA2058.csv")
a375 = read.csv("priorA375.csv")
skmel19 = read.csv("priorSKMEL19.csv")
skmel28 = read.csv("priorSKMEL28.csv")
skmel5 = read.csv("priorSubclones_Mix.csv")
wm88 = read.csv("priorWM88.csv")
wm164 = read.csv("priorWM164.csv")
wm793 = read.csv("priorWM793.csv")


# ===========================================================================================================
# 
# ===========================================================================================================
load("MCMC1_A375_150000_Run.RData") # A375
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data1 <- log(data[,c(1,2,3,4)])


# load("MCMC_Subclones_Mix_150000_Run.RData") # SKMEL5
# MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
# data <- data.frame(MCMC1$pars)
# data1 <- data[,c(1,2,3,4)]
# randomrows = sample(nrow(data1), 2000)
# store = data.frame()
# for(i in 1:length(randomrows)){
#   pars <- data1[i,]
#   dc = data.frame(k_rs = pars$k_ab, k_sr = pars$k_ba, k_es = pars$k_cb, k_se = pars$k_bc)
#   store = rbind(store, dc)
# }
# require(PerformanceAnalytics)
# chart.Correlation(store) 
# 
# library(psych)
# 
# pairs.panels(store, scale = F)
# # pairs(store)
# # pairs(store)
# panel.pearson <- function(x, y, ...) {
#   horizontal <- (par("usr")[1] + par("usr")[2]) / 2; 
#   vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
#   text(horizontal, vertical, format(abs(cor(x,y)), digits=2)) 
# }
# pairs(store[1:4], main = "Parameter Correlation", pch = 21, 
#       bg = c("red","green3","blue"), upper.panel=panel.pearson)
# ===========================================================================================================
load("MCMC1_SKMEL28_2e+05_Run.RData") # SKMEL28
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/1.6),]
data <- data.frame(MCMC1$pars)
data2 <- log(data[,c(1,2,3,4)])

load("MCMC1_SKMEL19_150000_Run.RData") # SKMEL19
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data3 <- log(data[,c(1,2,3,4)])


load("MCMC_Subclones_Mix_150000_Run.RData") # SKMEL5
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data4 <- log(data[,c(1,2,3,4)])


load("MCMC1_WM88_150000_Run.RData") # WM88
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data5 <- log(data[,c(1,2,3,4)])

load("MCMC1_WM793_150000_Run.RData") # WM793
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data6 <- log(data[,c(1,2,3,4)])

load("MCMC1_WM164_1e_05_Run.RData") # WM164
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data7 <- log(data[,c(1,2,3,4)])

load("MCMC1_A2058_150000_Run.RData") # A2058
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data8 <- log(data[,c(1,2,3,4)])


# ===========================================================================================================
# Plot the distributions of the rates on log scale
# ===========================================================================================================
par(ps = 12, cex = 1, cex.axis = 1)
# k_ab
pdf(paste0(output, "/", "A375_k_ab.pdf"), width=3.5, height=3.5)
hist(data1$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL28_k_ab.pdf"), width=3.5, height=3.5)
hist(data2$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL19_k_ab.pdf"), width=3.5, height=3.5)
hist(data3$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL5_k_ab.pdf"), width=3.5, height=3.5)
hist(data4$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM88_k_ab.pdf"), width=3.5, height=3.5)
hist(data5$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM793_k_ab.pdf"), width=3.5, height=3.5)
hist(data6$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM164_k_ab.pdf"), width=3.5, height=3.5)
hist(data7$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "A2058_k_ab.pdf"), width=3.5, height=3.5)
hist(data8$k_ab, breaks=seq(-16, -2, by=0.05), main="k_ab", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

# k_ba
pdf(paste0(output, "/", "A375_k_ba.pdf"), width=3.5, height=3.5)
hist(data1$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL28_k_ba.pdf"), width=3.5, height=3.5)
hist(data2$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL19_k_ba.pdf"), width=3.5, height=3.5)
hist(data3$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL5_k_ba.pdf"), width=3.5, height=3.5)
hist(data4$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM88_k_ba.pdf"), width=3.5, height=3.5)
hist(data5$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM793_k_ba.pdf"), width=3.5, height=3.5)
hist(data6$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "WM164_k_ba.pdf"), width=3.5, height=3.5)
hist(data7$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

pdf(paste0(output, "/", "A2058_k_ba.pdf"), width=3.5, height=3.5)
hist(data8$k_ba, breaks=seq(-16, -2, by=0.05), main="k_ba", freq = F, ylim=c(0,1.5), xlab="", xlim=c(-10,-2))
dev.off()

# k_cb
pdf(paste0(output, "/", "A375_k_cb.pdf"), width=3.5, height=3.5)
hist(data1$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL28_k_cb.pdf"), width=3.5, height=3.5)
hist(data2$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL19_k_cb.pdf"), width=3.5, height=3.5)
hist(data3$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL5_k_cb.pdf"), width=3.5, height=3.5)
hist(data4$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "WM88_k_cb.pdf"), width=3.5, height=3.5)
hist(data5$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "WM793_k_cb.pdf"), width=3.5, height=3.5)
hist(data6$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "WM164_k_cb.pdf"), width=3.5, height=3.5)
hist(data7$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-5,-2))
dev.off()

pdf(paste0(output, "/", "A2058_k_cb.pdf"), width=3.5, height=3.5)
hist(data8$k_cb, breaks=seq(-25, -2, by=0.05), main="k_cb", freq = F, ylim=c(0,5), xlab="", xlim=c(-15,0))
dev.off()

# k_bc
pdf(paste0(output, "/", "A375_k_bc.pdf"), width=3.5, height=3.5)
hist(data1$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL28_k_bc.pdf"), width=3.5, height=3.5)
hist(data2$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL19_k_bc.pdf"), width=3.5, height=3.5)
hist(data3$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "SKMEL5_k_bc.pdf"), width=3.5, height=3.5)
hist(data4$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "WM88_k_bc.pdf"), width=3.5, height=3.5)
hist(data5$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "WM793_k_bc.pdf"), width=3.5, height=3.5)
hist(data6$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "WM164_k_bc.pdf"), width=3.5, height=3.5)
hist(data7$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()

pdf(paste0(output, "/", "A2058_k_bc.pdf"), width=3.5, height=3.5)
hist(data8$k_bc, breaks=seq(-20, -2, by=0.05), main="k_bc", freq = F, ylim=c(0,2), xlab="", xlim=c(-15,-2))
dev.off()


# ===========================================================================================================
# ===========================================================================================================
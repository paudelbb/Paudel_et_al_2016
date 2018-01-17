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

cell = "skmel5"
output = paste0(output, "/", fig, "/", cell)
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)

files1 = grep("prior", dir(), value = T)

skmel5 = read.csv("priorSubclones_Mix.csv")

# ===========================================================================================================
# ===========================================================================================================

dt = skmel5
p1 = subset(dt, dt$k_ab >= 0 & dt$k_ba >= 0 & dt$k_bc >= 0 & dt$k_cb >= 0)
#p1 = subset(p1, p1$k_ab <= 0.06 & p1$k_ba <= 0.06 & p1$k_bc <= 0.06 & p1$k_cb <= 0.06)
a = subset(p1, p1$k_ab <= 0.06)
b = subset(p1, p1$k_ba <= 0.06)
c = subset(p1, p1$k_cb <= 0.06)
d = subset(p1, p1$k_bc <= 0.06)


# ===========================================================================================================
# MCMC runs
# ===========================================================================================================
load("MCMC_Subclones_Mix_150000_Run.RData") # SKMEL5
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
data <- data.frame(MCMC1$pars)
data1 <- log(data[,c(1,2,3,4)])
# ===========================================================================================================
x = 8
pdf(paste0(output, "/", cell, "k_rs.pdf"), width=3, height=3)
par(ps = x, cex = 1, cex.axis = 1)
hist(data1$k_ab,  freq = F, breaks=seq(-16, 0, by=0.25), col = "brown", xlim=c(-16, 0), ylim = c(0,1), main = "k_rs")
hist(log(a$k_ab), freq = F, breaks=seq(-16, 0, by=0.25), col = rgb(0,0,1,0.4),  add=T)
dev.off()

pdf(paste0(output, "/", cell, "k_sr.pdf"), width=3, height=3)
par(ps = x, cex = 1, cex.axis = 1)
hist(data1$k_ba,  freq = F, breaks=seq(-16, 0, by=0.25), col = "brown", xlim=c(-16, 0), ylim = c(0,1), main = "k_sr")
hist(log(b$k_ba), freq = F, breaks=seq(-16, 0, by=0.25), col = rgb(0,0,1,0.4),  add=T)
dev.off()

pdf(paste0(output, "/", cell, "k_es.pdf"), width=3, height=3)
par(ps = x, cex = 1, cex.axis = 1)
hist(data1$k_cb,  freq = F, breaks=seq(-16, 0, by=0.25), col = "brown", xlim=c(-16, 0), ylim = c(0,1), main = "k_es")
hist(log(c$k_cb), freq = F, breaks=seq(-16, 0, by=0.25), col = rgb(0,0,1,0.4),  add=T)
dev.off()

pdf(paste0(output, "/", cell, "k_se.pdf"), width=3, height=3)
par(ps = x, cex = 1, cex.axis = 1)
hist(data1$k_bc,  freq = F, breaks=seq(-20, 0, by=0.25), col = "brown", xlim=c(-16, 0), ylim = c(0,1), main = "k_se")
hist(log(d$k_bc), freq = F, breaks=seq(-20, 0, by=0.25), col = rgb(0,0,1,0.4),  add=T)
dev.off()
# ===========================================================================================================

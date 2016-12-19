# ===========================================================================================================
# Supp. Figure 6A: MCMC Traceplot for different cell lines
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
fig = "Supp.Fig6" # Which figure is this in paper?
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
pdf(paste0(output, "/", "A375_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()

load("MCMC_Subclones_Mix_150000_Run.RData") # SKMEL5
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "SKMEL5_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()


load("MCMC1_SKMEL19_150000_Run.RData") # SKMEL19
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "SKMEL19_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()

load("MCMC1_SKMEL28_150000_Run.RData") # SKMEL28
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "SKMEL28_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()

load("MCMC1_WM88_150000_Run.RData") # WM88
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "WM88_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()

load("MCMC1_WM164_1e_05_Run.RData") # WM164
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "WM164_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()

load("MCMC1_WM793_150000_Run.RData") # WM793
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
pdf(paste0(output, "/", "WM793_traceplot.pdf"), width=5.5, height=4)
plot(MCMC1)
dev.off()
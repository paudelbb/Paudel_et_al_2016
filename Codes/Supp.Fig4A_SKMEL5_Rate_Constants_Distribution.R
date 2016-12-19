# ===========================================================================================================
# Figure 4A: SKMEL5 distribution of rate constants
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
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)
chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4)])


# ===========================================================================================================
# Plot the rate constants distributions
# ===========================================================================================================
group = names(rates)
group <- noquote(group)
labels <- group
x = c(1,2,3,4)

pdf(paste0(output, "/Boxplot_of_SKMEL5_MCMC_rate_constants.pdf"), width=3.2, height=3.2)
boxplot(rates, notch=T, col=grey.colors(4), xaxt="n")
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
dev.off()
# ===========================================================================================================

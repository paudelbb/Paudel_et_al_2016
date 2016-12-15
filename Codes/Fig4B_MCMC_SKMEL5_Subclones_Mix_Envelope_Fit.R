# ===========================================================================================================
# Figure 4: SKMEL5 Subclones Mix Model Fitting and MCMC Algorithms
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
fig = "Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
#                                         Read experimental data
# ===========================================================================================================
data <- read.csv("SKMEL5_Subclones_Mix.csv", header=T, sep=",")
s1 <- data
d_exp <- data.matrix(s1[,c("Time", "nl2")]);
colnames(d_exp) <- c("time", "nl2")
d_exp1 <- data.frame(time=seq(0,300), nl2=c(-1000))
cell <- "Subclones_Mix"
# ===========================================================================================================
#                               SKMEL5 Subclones Mix MCMC Envelope
# ===========================================================================================================

dn <- read.csv(paste0(getwd(), "/MCMC_ensembles_", cell,".csv"))

# ===========================================================================================================
#               Plot experimental data with 1000 mcMC emsembles
# ===========================================================================================================
ggplot(data = data, aes(x=data$Time, y=data$nl2, col=CellLine)) +
  theme_bw() + geom_smooth(span=.3, aes(fill=CellLine), fill = "grey", alpha=0.6, data=data, 
  method = "loess", size=.5, level=0.95)+ 
  scale_colour_manual(values=c("black")) + ylim(-1, 1.0) + xlim(0,260)+
  theme(legend.position="none") + 
  geom_line(aes(x=dn$time, y=dn$nl2), data=dn, col="blue", alpha=0.5)+ ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=10)) + 
  ggsave(paste0(cell, " + ", "MCMC_envelope", ".png"), path=output, width=3, height=3)
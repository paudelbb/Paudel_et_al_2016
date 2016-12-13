#=====================================================================================================
# Figure 1: Population Level Response for different cell lines
#=====================================================================================================
#         Population level response of BRAF-mutated melanoma cell lines to PLX4720
#     Fraction of total cells that are positive for FUCCI--live cell reporter for cycling cells
#                   CellLines used: SKMEL5, A375, SKMEL28, WM164
#                 Concentrations for idling phenotype in cell lines:
#                 SKMEL5(8uM), A375(8uM), SKMEL28(16uM), WM164(32uM)
#=====================================================================================================

library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read the experimental data
# Define the time frame for the idling phenotype
# SKMEL5 (post 168hrs), A375 (post 180 hrs), WM164 (post 200hrs), SKMEL28 (168hrs)
#=====================================================================================================

data <- read.csv("Fucci_positive_percentage.csv", header=T, sep=",")

#=====================================================================================================
# SKMEL5 Percentage of Fucci positive cells
#=====================================================================================================

cell <- "SKMEL5"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)
s1 <- subset(s1, s1$Time>=168)

ggplot(data = s1, aes(x=s1$Time, y=s1$ffrac, col=Date))+ 
  theme_bw()+geom_smooth(span=.4, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(0,20)+ xlim(160, 350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") + 
  ggsave(paste0(cell, " + ", cnc, "μΜ_","Fucci_positive.pdf"), path=output, width=3, height=3)

#=====================================================================================================
# A375 Percentage of Fucci positive cells
#=====================================================================================================
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures")
fig = "Supp.Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
cell <- "A375"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)
s1 <- subset(s1, s1$Time>=170)

ggplot(data = s1, aes(x=s1$Time, y=s1$ffrac, col=Date))+ 
  theme_bw()+geom_smooth(span=.4, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(0,20)+ xlim(160, 350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") + 
  ggsave(paste0(cell, " + ", cnc, "μΜ_","Fucci_positive.pdf"), path=output, width=3, height=3)

#=====================================================================================================
# WM164 Percentage of Fucci positive cells
#=====================================================================================================
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures")
fig = "Supp.Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
cell <- "WM164"
cnc <- 32
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)
s1 <- subset(s1, s1$Time>=170)

ggplot(data = s1, aes(x=s1$Time, y=s1$ffrac, col=Date))+ 
  theme_bw()+geom_smooth(span=.4, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(0,20)+ xlim(160, 350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") + 
  ggsave(paste0(cell, " + ", cnc, "μΜ_","Fucci_positive.pdf"), path=output, width=3, height=3)


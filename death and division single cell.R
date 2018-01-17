#=====================================================================================================
# Figure 3B: Single cell response for SKMEL5 Single-cell-derived sublines in idling phase
#=====================================================================================================
#         Population level response of Sublines to PLX4720
#         CellLines used: SC01, SC07 and SC10
#         Concentrations used: 8uM
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read experimental data
#=====================================================================================================
data <- read.csv("SKMEL5_Subclones_SingleCell_Response.csv", header=T, sep=",")
data <- subset(data, data$phase=="idling")
div <- subset(data, data$Died==FALSE)
deg <- subset(data, data$Died==TRUE)
# boxplot of the time spent in G1 vs S-G2-M-time by subclones
boxplot(G1_time ~ CellLine, data=div, notch=T)
boxplot(G1_time ~ CellLine, data=deg)

boxplot(S_G2_M_time ~ CellLine, data=div, notch=T)
boxplot(S_G2_M_time ~ CellLine, data=deg, outline=F)

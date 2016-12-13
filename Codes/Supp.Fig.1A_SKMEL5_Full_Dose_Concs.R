#=====================================================================================================
# Figure 1: Population Level Response for different cell lines
#=====================================================================================================
#         Population level response of BRAF-mutated melanoma cell lines to PLX4720
#         CellLines used: SKMEL5, SKMEL19, SKMEL28, WM88, WM164, WM793, A375
#                 Concentrations for idling phenotype in cell lines:
#         SKMEL5(8uM), SKMEL19(8uM), SKMEL28(16uM), WM88(8uM), WM164(32uM), WM793(16uM), A375(8uM)
#=====================================================================================================
#=====================================================================================================
library("Hmisc")
library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Supp.Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
cell <- "SKMEL5"
s1 <- subset(data, data$CellLine==cell & data$Date =="20150617")
s1$Date <- as.character(s1$Date)

#=====================================================================================================
#  Summarize the data and get mean+-sd
#=====================================================================================================

source("~/Paudel_et_al_2016/Codes/summarySE.R")
d <- s1[s1$conc!=0.25, ]
d <- d[d$conc!=0.125, ]
df <- summarySE(d, measurevar="nl2", groupvars=c("Time", "conc"))
d <- df

#=====================================================================================================
# Graph function
#=====================================================================================================

par(ps = 12, cex = 1, cex.axis = 1)
graph <- function(time, nl2, conc, title, sub)
{
  y_min = -1
  y_max = 5.5
  plot(  time, nl2, 
         xlim=c(0,(max(d$Time)+35)), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub)
  with (
    data = d, expr = errbar(Time, nl2, nl2+se, nl2-se, add=T, pch=21, cap=.002))
  for(tc in unique(conc))
  {
    ccc <- conc[tc==conc]
    lines(time[tc==conc], nl2[tc==conc])
    text((max(d$Time)-5), nl2[tc==conc][length(ccc)], label=paste0(tc,'µM'), pos=4, cex=1.0)
  }
}

#=====================================================================================================
# Plot the full dose range in SKMEL5
#=====================================================================================================

pdf(paste0(output, "/SKMEL5_full_PLX4720_Range.pdf"), height=6, width=4)
graph(d$Time, d$nl2, d$conc, title="", sub="")
dev.off()


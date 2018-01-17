##################################################################################
# SKMEL28 Parental + PLX4720 response
##################################################################################
dd <- "X:/Subclones_Bishal/LongTerm"
setwd("X:/Subclones_Bishal/LongTerm")

d <- read.csv("20160203 WM cells Processed.csv", header=T, sep=",")
##################################################################################
library(ggplot2)
source("X:\\Users\\paudelbb\\Paper_1\\Used R Codes\\summarySE.R")
df <- summarySE(d, measurevar="nl2", groupvars=c("Time","drug", "conc", "CellLine"))
d <- df
##################################################################################
##################################################################################
library("Hmisc")
graph <- function(data, time, nl2, conc, title, sub)
{
  y_min = -1
  y_max = 6
  d = data
  plot(  time, nl2, 
         xlim=c(0,(max(d$Time)+20)), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub)
  with (
    data = d, expr = errbar(Time, nl2, nl2+sd, nl2-sd, add=T, pch=16, cap=.002))
  for(tc in unique(conc))
  {
    ccc <- conc[tc==conc]
    lines(time[tc==conc], nl2[tc==conc])
    text((max(d$Time)-2), nl2[tc==conc][length(ccc)], label=paste0(tc,'uM'), pos=4, cex=1.0)
  }
}
d1 <- subset(d, d$CellLine=="WM164")
d2 <- subset(d, d$CellLine=="WM793")
d3 <- subset(d, d$CellLine=="WM88")

par(mfrow=c(1,3))
graph(d1, d1$Time, d1$nl2, d1$conc, title="", sub="")
graph(d2, d2$Time, d2$nl2, d2$conc, title="", sub="")
graph(d3, d3$Time, d3$nl2, d3$conc, title="", sub="")

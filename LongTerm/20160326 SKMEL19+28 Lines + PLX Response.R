#=======================================================================================
#
# Analysis of Subclones response -- Untreated vs Post idle in different inhibitors
# Compares untreated vs postidle subclones--three different concentrations and a control
#
#=======================================================================================
# Set the working directory, read the file, and pre-process the data
require(ggplot2)
library("Hmisc")
source("X:\\Users\\paudelbb\\Paper_1\\Used R Codes\\summarySE.R")
dd <- "X:/Data/CELLAVISTA/Bishal Paudel/Subclones_Bishal/LongTerm"
setwd(dd)
output <- "X:/Data/CELLAVISTA/Bishal Paudel/Subclones_Bishal/LongTerm"
setwd(dd)
d1 <- read.csv("20160326 SKMEL19+28 + PLX4720 Processed.csv", header=T, sep=",")
d2 <- read.csv("20150227 SKMEL28 Processed.csv", header=T, sep=",")
# 
# #=======================================================================================
# a <- "SKMEL28"
# cnc <- 16
# s1 <- subset(d1, d1$CellLine==a & d1$conc==cnc)
# s2 <- subset(d2, d2$CellLine==a & d2$conc==cnc)
# s1 <- rbind(s1, s2)
# s1$Date <- as.character(s1$Date)
# ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
#   theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
#   scale_colour_manual(values=c("blue")) + ylim(-3.0,3.0)+ xlim(0,500)+
#   theme(legend.position="none") + theme(axis.text=element_text(size=16))+theme(text = element_text(size=16))+
#   ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
# + ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)
# 


a <- "SKMEL19"
cnc <- 4
s1 <- subset(d, d$CellLine==a & d$conc==cnc)
ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=CellLine))+
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
  scale_colour_manual(values=c("blue")) + ylim(-3.0,3.0)+ xlim(-0,350)+
  theme(legend.position="none")+ theme(axis.text=element_text(size=16))+theme(text = element_text(size=16))+
  ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
+ ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)

















# 
# 
# 
# df <- summarySE(d, measurevar="nl2", groupvars=c("Time", "conc", "CellLine"))
# d <- df

# #=======================================================================================
# graph <- function(data, time, nl2, conc, title, sub)
# {
#   d <- data
#   y_min = -3
#   y_max = 6
#   plot(  time, nl2, 
#          xlim=c(0,(max(d$Time)+60)), ylim=c(y_min, y_max),
#          xlab='',ylab='',
#          main=title, sub=sub, cex.axis=1.5)
#   with (
#     data = d, expr = errbar(Time, nl2, nl2+se, nl2-se, add=T, pch=16, cap=.002))
#   for(tc in unique(conc))
#   {
#     ccc <- conc[tc==conc]
#     lines(time[tc==conc], nl2[tc==conc])
#     text((max(d$Time)-2), nl2[tc==conc][length(ccc)], label=paste0	(''), pos=4)
#   }
# }
# #=======================================================================================
# #d <- subset(d, d$conc==16)
# d1 <- subset(d, d$CellLine==unique(d$CellLine)[1])
# d2 <- subset(d, d$CellLine==unique(d$CellLine)[2])
# 
# par(mfrow=c(1,2))
# graph(d1, d1$Time, d1$nl2, d1$conc, title="", sub="")
# graph(d2, d2$Time, d2$nl2, d2$conc, title="", sub="")
# 

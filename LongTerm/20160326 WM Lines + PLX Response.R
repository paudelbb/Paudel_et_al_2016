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
dd <- "X:/Subclones_Bishal/LongTerm"
output <- "X:/Subclones_Bishal/LongTerm"
setwd(dd)
d <- read.csv("20160326 WM Lines + PLX4720 Processed.csv", header=T, sep=",")

#=======================================================================================
a <- "WM88"
cnc <- 8
s1 <- subset(d, d$CellLine==a & d$conc==cnc)
 ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=CellLine))+
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
  scale_colour_manual(values=c("red")) + ylim(-4.5,5.0)+ xlim(-0,350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=20))+
  ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
#+ ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)

cnc <- 16
s2 <- subset(d, d$CellLine==a & d$conc==cnc)
plot2 <- ggplot(data = s2, aes(x=s2$Time, y=s2$nl2, col=CellLine)) +geom_point(aes(shape=CellLine))+
  theme_bw()+geom_smooth(span=.3, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
  scale_colour_manual(values=c("red")) + ylim(-4.5,5.5)+ xlim(-0,350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=20))+
  ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
#+ ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)

cnc <- 8
s3 <- subset(d, d$CellLine==a & d$conc==cnc)
plot3 <- ggplot(data = s3, aes(x=s3$Time, y=s3$nl2, col=CellLine)) +geom_point(aes(shape=CellLine))+
  theme_bw()+geom_smooth(span=.3, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
  scale_colour_manual(values=c("red")) + ylim(-4.5,5.5)+ xlim(-0,350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=20))+
  ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
#+ ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)


cnc <- 2
s4 <- subset(d, d$CellLine==a & d$conc==cnc)
plot4 <- ggplot(data = s4, aes(x=s4$Time, y=s4$nl2, col=CellLine)) +geom_point(aes(shape=CellLine))+
  theme_bw()+geom_smooth(span=.3, aes(group=1), method = "loess", size=1, alpha=0.7)+ 
  scale_colour_manual(values=c("red")) + ylim(-4.5,5.5)+ xlim(-0,350)+
  theme(legend.position="none") + theme(axis.text=element_text(size=20))+
  ggtitle(paste0(a, " + ", cnc, "uM"))+ labs(x="", y="")
#+ ggsave(paste0(a, " + ", cnc, "uM", ".pdf"), path=output, width=4.5, height=4.5)
# grid.arrange(plot4,plot3, plot2, plot1,  ncol=2, nrow=2)












# 
# 
# 
# 
# df <- summarySE(d, measurevar="nl2", groupvars=c("Time", "conc", "CellLine"))
# d <- df
# 
# #=======================================================================================
# graph <- function(data, time, nl2, conc, title, sub)
# {
#   d <- data
#   y_min = -1
#   y_max = 3
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
# d1 <- subset(d, d$CellLine=="WM164")
# d2 <- subset(d, d$CellLine=="WM793")
# d3 <- subset(d, d$CellLine=="WM88")
# #par(mfrow=c(1,3))
# graph(d1, d1$Time, d1$nl2, d1$conc, title="", sub="")
# graph(d2, d2$Time, d2$nl2, d2$conc, title="", sub="")
# graph(d3, d3$Time, d3$nl2, d3$conc, title="", sub="")
# 
# 
# # abline(h=2.5, col="red")
# # abline(h=3.1, col="red")
# # abline(h=0, col="red")
# # abline(h=1.6, col="red")
# 
# # s1 <- subset(d2, d2$Time<150)
# # s2 <- subset(d2, d2$Time>150)
# # 
# # 
# # coef(lm(nl2~Time, data=s2))

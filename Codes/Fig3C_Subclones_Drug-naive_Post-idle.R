#=====================================================================================================
# Figure 3C: SKMEL5 Subclones response--drug-naive vs post-idle in different inhibitors
# Inhibitors used: Trametinib, Oligomycin, Cisplatin, BKM120, BEZ235
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)

setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)
ps = 12

source("~//Paudel_et_al_2016//Codes//summarySE.R") # Summarize function to get the summary statistics
source("~//Paudel_et_al_2016//Codes//graphFunction1.R") # main graph function (graph1)
source("~//Paudel_et_al_2016//Codes//graphFunction2.R") # graph function to add to an existing plot (graph2)

#=====================================================================================================
# Read experimental data
#=====================================================================================================

data <- read.csv("SKMEL5 Subclones Unt vs Postidle inhibitors.csv", header=T, sep=",")

s1 <- subset(data,  data$Date =="20160315" & data$Time<=100)
s1$Date <- as.character(s1$Date)

#=====================================================================================================
# Summarize data by time, concentration, cellline, conditions--mean+-sd_or_se
#=====================================================================================================

dt <- summarySE(s1, measurevar = "nl2", groupvars = c("CellLine", "Time", "drug", "cond", "conc"))

#=====================================================================================================
# Trametinib Response in drug-naive vs post-idle cells
#=====================================================================================================
inh <- "Trametinib"
dt1 <- subset(dt, dt$drug==inh)
d1 <- subset(dt1, dt1$cond=="Untreated" & dt1$conc==1)
d2 <- subset(dt1, dt1$cond=="Postidle"  & dt1$conc==1)
pdf(paste0(output, "/SKMEL5_Subclones_Trametinib.pdf"), width=3, height=4.2)
par(ps = ps, cex = 1, cex.axis = 1)
graph1(d1, d1$Time, d1$nl2, d1$CellLine, title="", sub="", col=c("black"), 0, 120, -2.5, 1)
graph2(d2, d2$Time, d2$nl2, d2$CellLine, title="", sub="", col=c("red"), 0, 120, -2.5, 1)
dev.off()
#=====================================================================================================
# Cisplatin Response in drug-naive vs post-idle cells
#=====================================================================================================
inh <- "Cisplatin"
dt1 <- subset(dt, dt$drug==inh)
d1 <- subset(dt1, dt1$cond=="Untreated" & dt1$conc==5)
d2 <- subset(dt1, dt1$cond=="Postidle"  & dt1$conc==5)
pdf(paste0(output, "/SKMEL5_Subclones_Cisplatin_early.pdf"), width=3, height=4.2)
par(ps = ps, cex = 1, cex.axis = 1)
graph1(d1, d1$Time, d1$nl2, d1$CellLine, title="", sub="", col=c("black"), 0, 120, -2.5, 1)
graph2(d2, d2$Time, d2$nl2, d2$CellLine, title="", sub="", col=c("red"), 0, 120, -2.5, 1)
dev.off()
#=====================================================================================================
# BKM120 response in drug-naive vs post-idle cells
#=====================================================================================================
data <- read.csv("SKMEL5 Subclones Unt vs Postidle inhibitors.csv", header=T, sep=",")
s1 <- subset(data,  data$Date =="20160811" & data$Time<=100)
s1$Date <- as.character(s1$Date)
dt <- summarySE(s1, measurevar = "nl2", groupvars = c("CellLine", "Time", "drug", "cond", "conc"))
inh <- "BKM120"
dt1 <- subset(dt, dt$drug==inh)
d1 <- subset(dt1, dt1$cond=="Untreated" & dt1$conc==1)
d2 <- subset(dt1, dt1$cond=="Postidle"  & dt1$conc==1)
pdf(paste0(output, "/SKMEL5_Subclones_BKM120.pdf"), width=3, height=4.2)
par(ps = ps, cex = 1, cex.axis = 1)
graph1(d1, d1$Time, d1$nl2, d1$CellLine, title="", sub="", col=c("black"), 0, 65, -2.5, 1)
graph2(d2, d2$Time, d2$nl2, d2$CellLine, title="", sub="", col=c("red"), 0, 65, -2.5, 1)
dev.off()

#=====================================================================================================
# BEZ235 Response in drug-naive vs post-idle cells
#=====================================================================================================

output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Supp.Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)

data <- read.csv("SKMEL5 Subclones Unt vs Postidle inhibitors.csv", header=T, sep=",")
s1 <- subset(data,  data$Date =="20160315" & data$Time<=100)
s1$Date <- as.character(s1$Date)
dt <- summarySE(s1, measurevar = "nl2", groupvars = c("CellLine", "Time", "drug", "cond", "conc"))
inh <- "BEZ-235"
dt1 <- subset(dt, dt$drug==inh)
d1 <- subset(dt1, dt1$cond=="Untreated" & dt1$conc==2)
d2 <- subset(dt1, dt1$cond=="Postidle"  & dt1$conc==2)
pdf(paste0(output, "/SKMEL5_Subclones_BEZ235.pdf"), width=3, height=4.2)
par(ps = ps, cex = 1, cex.axis = 1)
graph1(d1, d1$Time, d1$nl2, d1$CellLine, title="", sub="", col=c("black"), 0, 120, -2.5, 1)
graph2(d2, d2$Time, d2$nl2, d2$CellLine, title="", sub="", col=c("red"), 0, 120, -2.5, 1)
dev.off()
#=====================================================================================================
#=====================================================================================================
pdf(paste0(output,"/","Empty_legend_inhibitors.pdf"), width=4, height=4)
par(ps = ps, cex = 1, cex.axis = 1)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(0.55,1, c("Drug-naive", "Post-idle"), col = c("black","red"),horiz=F,
       lty=c(1,2), lwd=1.5, bg = "gray90", bty="n")
dev.off()

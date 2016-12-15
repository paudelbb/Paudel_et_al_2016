#=====================================================================================================
# Supp. Figure 3D: Single cell response for SKMEL5 Single-cell-derived sublines in initial phase
# 0-120 hours of treatment with 8uM PLX4720
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
fig = "Supp.Fig3" # Which figure is this in paper?
output = paste0(output, "/", fig)

#=====================================================================================================
# Read experimental data
#=====================================================================================================
data <- read.csv("SKMEL5_Subclones_SingleCell_Response.csv", header=T, sep=",")
data <- subset(data, data$phase=="early")
data <- subset(data, data$Birth.time..h.<=120)
# deg is for death events, div is for division events, doe is for end of experiment events
deg <- subset(data, data$Died==TRUE)
div <- subset(data, data$Died==FALSE)
div <- subset(div, div$EOE!=TRUE)
doe <- subset(data, data$EOE==TRUE)

#=====================================================================================================
# Response for different sublines
#=====================================================================================================
set.seed(12345)
# death events in sublines
t1 = subset(deg, deg$CellLine=="SC01"); t2 = subset(deg, deg$CellLine=="SC07"); t3 = subset(deg, deg$CellLine=="SC10")
l1 = length(t1$X); l2 = length(t2$X); l3 = length(t3$X)
t_1 = t1
t_2 = t2
t_3 = t3
# division events in sublines
s1 = subset(div, div$CellLine=="SC01"); s2 = subset(div, div$CellLine=="SC07"); s3 = subset(div, div$CellLine=="SC10")
k1 = length(s1$X); k2 = length(s2$X); k3 = length(s3$X)
s_1 = s1
s_2 = s2
s_3 = s3
# EOE events in sublines
r1 = subset(doe, doe$CellLine=="SC01"); r2 = subset(doe, doe$CellLine=="SC07"); r3 = subset(doe, doe$CellLine=="SC10")
m1 = length(r1$X); m2 = length(r2$X); m3 = length(r3$X)
r_1 = r1
r_2 = r2
r_3 = r3

deg <- rbind(t_1, t_2, t_3)
div <- rbind(s_1, s_2, s_3)
doe <- rbind(r_1, r_2, r_3)
#=====================================================================================================
#         Plot single cell responses
#=====================================================================================================

wd = 3; ht =3 
# size of the texts on axis and on plots
ps = 8

pdf(paste0(output,"/","SKMEL5 Subclone01 + SCT_Initial.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC01")),], 
     col='red', pch=4, cex=1, xlim=c(0,120), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC01")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC01")),], col='black', pch=22, cex=1.0)
dev.off()

pdf(paste0(output,"/","SKMEL5 Subclone07 + SCT_Initial.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC07")),], 
     col='red', pch=4, cex=1, xlim=c(0,120), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC07")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC07")),], col='black', pch=22, cex=1.0)
dev.off()

pdf(paste0(output,"/","SKMEL5 Subclone10 + SCT_Initial.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC10")),], 
     col='red', pch=4, cex=1, xlim=c(0,120), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC10")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC10")),], col='black', pch=22, cex=1.0)
dev.off()

#=====================================================================================================
#=====================================================================================================


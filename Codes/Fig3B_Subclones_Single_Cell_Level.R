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
t_1 = data.frame(t1[sample(nrow(t1), min(l1, l2, l3), replace=F), ])
t_2 = data.frame(t2[sample(nrow(t2), min(l1, l2, l3), replace=F), ])
t_3 = data.frame(t3[sample(nrow(t3), min(l1, l2, l3), replace=F), ])
# division events in sublines
s1 = subset(div, div$CellLine=="SC01"); s2 = subset(div, div$CellLine=="SC07"); s3 = subset(div, div$CellLine=="SC10")
k1 = length(s1$X); k2 = length(s2$X); k3 = length(s3$X)
s_1 = data.frame(s1[sample(nrow(s1), min(k1, k2, k3), replace=F), ])
s_2 = data.frame(s2[sample(nrow(s2), min(k1, k2, k3), replace=F), ])
s_3 = data.frame(s3[sample(nrow(s3), min(k1, k2, k3), replace=F), ])
# EOE events in sublines
r1 = subset(doe, doe$CellLine=="SC01"); r2 = subset(doe, doe$CellLine=="SC07"); r3 = subset(doe, doe$CellLine=="SC10")
m1 = length(r1$X); m2 = length(r2$X); m3 = length(r3$X)
r_1 = data.frame(r1[sample(nrow(r1), min(m1, m2, m3), replace=F), ])
r_2 = data.frame(r2[sample(nrow(r2), min(m1, m2, m3), replace=F), ])
r_3 = data.frame(r3[sample(nrow(r3), min(m1, m2, m3), replace=F), ])

deg <- rbind(t_1, t_2, t_3)
div <- rbind(s_1, s_2, s_3)
doe <- rbind(r_1, r_2, r_3)
#=====================================================================================================
#         Plot single cell responses
#=====================================================================================================

wd = 3; ht =3 
ps = 8

pdf(paste0(output,"/","SKMEL5 Subclone01 + SCT_Idling.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC01")),], 
     col='red', pch=4, cex=1, xlim=c(160,260), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC01")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC01")),], col='black', pch=22, cex=1.0)
dev.off()

pdf(paste0(output,"/","SKMEL5 Subclone07 + SCT_Idling.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC07")),], 
     col='red', pch=4, cex=1, xlim=c(160,260), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC07")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC07")),], col='black', pch=22, cex=1.0)
dev.off()

pdf(paste0(output,"/","SKMEL5 Subclone10 + SCT_Idling.pdf"), width=wd, height=ht)
par(ps = ps, cex = 1, cex.axis = 1)
plot(Lifetime..h.~Birth.time..h., data=deg[(deg$CellLine %in% c("SC10")),], 
     col='red', pch=4, cex=1, xlim=c(160,260), ylim=c(0,120), xlab="", ylab="")
points(Lifetime..h.~Birth.time..h., data=div[(div$CellLine %in% c("SC10")),], pch=18, cex=1.0)
points(Lifetime..h.~Birth.time..h., data=doe[(doe$CellLine %in% c("SC10")),], col='black', pch=22, cex=1.0)
dev.off()

#=====================================================================================================
#=====================================================================================================
pdf(paste0(output,"/","Empty_legend.pdf"), width=4, height=4)
par(ps = ps, cex = 1, cex.axis = 1)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(0.55,1, c("Division", "Death", "EOE"), col = c("black","red","black"),horiz=F,
       pch = c(18,4, 22),bg = "gray90", bty="n")
dev.off()

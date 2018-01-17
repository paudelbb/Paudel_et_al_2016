#==============================================================================================
# SKMEL5 response in vemurafenib and dabrafenib
#==============================================================================================
dd <- "/Users/paudelbb/Dropbox/PhD Work/Melanoma 2017/Vem_Dab"
setwd(dd)
output <- paste0(dd, "/Figs")
#==============================================================================================
# Read the dataset from the folder
data1 <- read.csv("20170301 SKMEL5 parental vemurafenib and dabrafenib processed_1_2.csv")
data2 <- read.csv("20170312 SKMEL5 parental vemurafenib and dabrafenib processed_2_2.csv")
data <- rbind(data1, data2)
data <- subset(data, data$Time <= 380)
#==============================================================================================
# Which drug and what concentration to pick?
dr <- "dabrafenib"; 
if(dr=="vemurafenib") cnc = 16 else if(dr=="dabrafenib") cnc=4
s1 <- subset(data, data$conc==cnc & data$drug==dr)
s1$Date <- as.character(s1$Date)
s1$conc <- as.character(s1$conc)

#==============================================================================================
# Plot the drug response for the drug and concentration desired;
#==============================================================================================
ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=conc)) +
  theme_bw()+ 
  geom_smooth(span=.25, aes(group=1),data=s1, method = "loess", size=.5, alpha=0.8, col="blue")+ 
  scale_colour_manual(values=rainbow(1)) + ylim(-2, 6)+
  theme(legend.position="none") +
  ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12)) +
  ggsave(paste0("SKMEL5 Parental_", dr, "_", cnc, "uM.pdf"), path=output, width=3, height=3)

#==============================================================================================
#==============================================================================================
# How do the rates of proliferation change in early vs late phase in drug
#==============================================================================================
lmd <- function(l2){
  timePt <- unique(s$Time)
  coef(lm(l2~timePt))[2]
}
#=============================================================================================
# Calcuate the DIP rates 
#==================Early Phase===============================================
if(dr=="vemurafenib") cnc = 16 else if(dr=="dabrafenib") cnc=4
if(dr=="vemurafenib") time1 = 60 else if(dr=="dabrafenib") time1 = 20
s1 <- subset(data, data$conc==cnc & data$drug==dr)
s <- subset(s1,  s1$Time>=time1 & s1$Time <= 144)
dn <- data.frame()
for(i in unique(s$Date)){
  temp <- subset(s, s$Date==i)
  for(j in unique(temp$rep)){
    temp1 <- subset(temp, temp$rep==j)
    dc <- data.frame(Date = unique(temp1$Date), 
                     CellLine = unique(temp1$CellLine), drug = unique(temp1$drug), conc=unique(temp1$conc),
                     rep = unique(temp1$rep), rates=coef(lm(nl2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
  }
}
rate1 <- dn
#==================Late Phase===============================================
s <- subset(s1,  s1$Time>=170)
dn <- data.frame()
for(i in unique(s$Date)){
  temp <- subset(s, s$Date==i)
  for(j in unique(temp$rep)){
    temp1 <- subset(temp, temp$rep==j)
    dc <- data.frame(Date = unique(temp1$Date), 
                     CellLine = unique(temp1$CellLine), drug = unique(temp1$drug), conc=unique(temp1$conc),
                     rep = unique(temp1$rep), rates=coef(lm(nl2~Time, data=temp1))[2])
    dn <- rbind(dn, dc)
  }
}
rate2 <- dn
rate1$phase <- "early"
rate2$phase <- "late"
nrate <- rbind(rate1, rate2)
ttest = t.test(rate1$rates, rate2$rates)
statp <- round(ttest$p.value, 4)
#=============================================================================================
# Plot the drug-induced rates in early and late phase;
source("/Users/paudelbb/Paudel_et_al_2016/Codes/summarySE.R")
df <- summarySE(nrate, measurevar="rates", groupvars=c("phase"))
mybarcol <- "gray20"
pdf(paste0(output,"/drug-induced_rates_early_vs_late_", dr, "_", cnc, "uM.pdf"), width=4, height=4)
mp <- barplot2(df$rates, beside=TRUE,
               col=grey.colors(2), ylab="",
               main = "", font.main = 5,
               sub = "", col.sub = mybarcol, ylim=c(-0.00, 0.020),
               cex.names = 0.8, plot.ci = TRUE, ci.l = df$rates-df$se, ci.u = df$rates+df$se,
               plot.grid = FALSE, cex=0.8)
group <- c(as.character(unique(df$phase)))
group <- noquote(group)
labels <- paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
dev.off()

#=============================================================================================
# Plot the drug responses along with DMSO control
#=============================================================================================

data1 <- read.csv("20170301 SKMEL5 parental vemurafenib and dabrafenib processed_1_2.csv")
data2 <- read.csv("20170312 SKMEL5 parental vemurafenib and dabrafenib processed_2_2.csv")
data <- rbind(data1, data2)
data <- subset(data, data$Time <= 380)
s1 <- subset(data, data$drug==dr)
if(dr=="vemurafenib") cnc = c(0,16) else if(dr=="dabrafenib") cnc= c(0,4)
s1 <- s1[(s1$conc %in% cnc), ]
dmso <- subset(s1, s1$conc==cnc[1] & s1$Time <= 180)
drug <- subset(s1, s1$conc==cnc[2])
s1 <- rbind(dmso, drug)
s1$conc <- as.character(s1$conc)
#==============================================================================================
ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=conc))+ 
  theme_bw()+geom_smooth(span=.30, aes(fill=conc), method = "loess", size=.5, alpha=0.7, col="blue")+ 
  scale_colour_manual(values=c("red", "blue")) + ylim(-2,6)+ xlim(0,370)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="")+
  ggsave(paste0("SKMEL5 Parental_", dr, "_", cnc[1],"_vs_",cnc[2], "uM.pdf"), path=output, width=3, height=3)
  
#==============================================================================================
# Plot the number of cells in DMSO control vs drug treatment
#==============================================================================================
data1 <- read.csv("20170301 SKMEL5 parental vemurafenib and dabrafenib processed_1_2.csv")
data2 <- read.csv("20170312 SKMEL5 parental vemurafenib and dabrafenib processed_2_2.csv")
data <- rbind(data1, data2)
data <- subset(data, data$Time <= 380)
# Which drug and which concentration to pick;
if(dr=="vemurafenib") cnc = c(0,16) else if(dr=="dabrafenib") cnc= c(0,4)
s1 <- s1[(s1$conc %in% cnc), ]
s1 <- subset(data, data$drug==dr)
s1 <- s1[(s1$conc %in% cnc), ]
idling <- c(207.7864, 206.7036)
s1 <- s1[(s1$Time %in% idling),]
#==============================================================================================
cellN <- summarySE(s1, measurevar="Cell.Nucleus", groupvars=c("conc"))

pdf(paste0(output,"/Cell_count_in_idling_", dr, "_", cnc[1], "_vs_", cnc[2], "uM.pdf"), width=4, height=4)
mp <- barplot2(cellN$Cell.Nucleus, beside=TRUE,
               col=grey.colors(2), ylab="",
               main = "", font.main = 5,
               sub = "", col.sub = mybarcol, ylim=c(0, 25000),
               cex.names = 0.8, plot.ci = TRUE, ci.l = cellN$Cell.Nucleus-cellN$se, ci.u = cellN$Cell.Nucleus+cellN$se,
               plot.grid = FALSE, cex=0.8)
group <- c(as.character(unique(cellN$conc)))
group <- noquote(group)
labels <- paste0(group, "uM")
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
dev.off()
ttest1 = t.test(s1$Cell.Nucleus[s1$conc==cnc[1]], s1$Cell.Nucleus[s1$conc==cnc[2]])
stat1 <- ttest1$p.value
#==============================================================================================

#=====================================================================================================
# Figure 2: Population Level Response for SKMEL5 Single-cell-derived sublines
#=====================================================================================================
# Population level response of Sublines to PLX4720, CellLines used: SC01, SC07 and SC10, conc: 8uM
#=====================================================================================================
#=====================================================================================================
library(gplots)
require(ggplot2)
require(Hmisc)
source("~//Paudel et al. 2017 Paper 1//SourceCode//summarySE.R") # Summarize function to get the summary statistics;
setwd("~/Paudel et al. 2017 Paper 1/SourceData") # Set the working directory to the folder with data;
# output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") # Set the output for plots
# fig = "Fig2" # Which figure is this in paper?
# output = paste0(output, "/", fig)

output <- "/Users/paudelbb/Dropbox/Paudel_et_al_2017_idle_paper" 
fig = "DefenseTalk" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("SKMEL5_Subclones_Response.csv", header=T, sep=",")


#=====================================================================================================
# SC01 pop Level Response
#=====================================================================================================
cell <- "SKMEL5 Subclone01"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.45, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)
#=====================================================================================================
cell <- "SKMEL5 Subclone07"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.45, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


cell <- "SKMEL5 Subclone10"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.45, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


#=====================================================================================================
# Plot the response
#=====================================================================================================
ggplot(data = data, aes(x=data$Time, y=data$nl2, col=CellLine)) +
  theme_bw()+ 
  geom_smooth(span=.5, aes(fill=CellLine),data=data, method = "loess", size=.5, alpha=0.6)+ 
  scale_colour_manual(values=c("red", "green", "deepskyblue")) + ylim(-1.5, 1.5)+ xlim(0, 350)+
  theme(legend.position="none") +
  ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12)) +
  ggsave("SKMEL5 Subclones + 8uM PLX4720.pdf", path=output, width=2.5, height=3)

#=====================================================================================================
#=====================================================================================================
# early
#=====================================================================================================
spt = 72
st = 15
dn = data.frame()
for(i in unique(data$Date)){
  temp = subset(data, data$Date==i)
  for(j in unique(temp$CellLine)){
    t = subset(temp, temp$CellLine==j)
    for(k in unique(t$conc)){
      t1 = subset(t, t$conc==k)
      for(l in unique(t1$Well)){
        t2 = subset(t1, t1$Well==l)
        early = subset(t2, t2$Time >=st & t2$Time <=spt)
        dc = data.frame(Date = unique(t2$Date), CellLine = unique(t2$CellLine), 
                        drug = unique(t2$drug), conc = unique(t2$conc), 
                        rate1 = coef(lm(nl2~Time, data=early))[2])
        dn = rbind(dn, dc)
      }
    }
  }
}
early = dn
#=====================================================================================================
# late
#=====================================================================================================
stemp = subset(data, data$Time >=100)
spt = 100
dn = data.frame()
for(i in unique(stemp$Date)){
  temp = subset(stemp, stemp$Date==i)
  for(j in unique(temp$CellLine)){
    t = subset(temp, temp$CellLine==j)
    for(k in unique(t$conc)){
      t1 = subset(t, t$conc==k)
      for(l in unique(t1$Well)){
        t2 = subset(t1, t1$Well==l)
        late = subset(t2, t2$Time >=spt)
        dc = data.frame(Date = unique(t2$Date), CellLine = unique(t2$CellLine), 
                        drug = unique(t2$drug), conc = unique(t2$conc), 
                        rate1 = coef(lm(nl2~Time, data=late))[2])
        dn = rbind(dn, dc)
      }
    }
  }
}
late = dn
#=====================================================================================================
#=====================================================================================================
early$phase = "initial"
late$phase  = "late"
dtn = rbind(early, late)
ds <- summarySE(dtn, measurevar="rate1", groupvars=c("CellLine", "conc", "phase"))


#pdf(paste0(output,"/","SKMEL5 Sublines Short-term DIP rates.pdf"), width=3, height=3)
par(ps = 10, cex = 1, cex.axis = 1)
mp = barplot2(ds$rate1, beside = T,
              col=grey.colors(nrow(ds)), ylab="DIP rates (doublings/h)",
              main="", font.main=5, sub = "", col.sub = mybarcol, ylim=c(-0.03, 0.02),
              plot.ci = T, ci.l = ds$rate1-ds$se, ci.u = ds$rate1+ds$se,
              plot.grid = F)
group = c(as.character((ds$CellLine)))
group = noquote(group)
labels = paste0("SC",group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1, 1.1), xpd = T, cex = 0.5)




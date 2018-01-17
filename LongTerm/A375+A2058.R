#=====================================================================================================================
# Analysis of A2058 response in PLX4720
#=====================================================================================================================
# Set the working directory, read the file, and pre-process the data
dd <- "/Users/paudelbb/Paudel_et_al_2016/LongTerm"
setwd(dd)
expName <- "20150227 A375+A2058+PLX4720.csv"
source("loadCellCount.R")
d <- loadCellCount(dd, 60)
#=======================================================================================
# Assign the cell type, different concentration, drugs etc 
cell <- c('A375','A2058')
d$CellLine <- rep(cell, each=3)
d$drug <- "plx"
d$conc <- rep(8)
#=======================================================================================
dn <- data.frame()
for(l in unique(d$CellLine)){
  s1 <- subset(d, d$CellLine==l)
  for(i in unique(s1$conc)){
    temp <- subset(s1, s1$conc==i)
    temp$rep <- rep(seq(1:(length(temp$Time)/length(unique(temp$Time)))))
    dn <- rbind(dn, temp)
  }
}
d <- dn
#===============================================================

#===============================================================
well <- as.character(unique(d$Well))
dn <- data.frame()
dn <- data.frame()
for (i in unique(d$Well)){
  temp <- subset(d, d$Well==i)
  temp$nl2 <- temp$l2 - temp$l2[temp$Time==0]
  dn <- rbind(dn, temp)
}
d <- dn
#=======================================================================================
d$Date <- "20150227"
d <- d[, c("Date","Row","Well","CellLine","Cell.Nucleus", "Time", "drug","conc", "rep", "l2", "nl2" )]
write.csv(d, file="20150227 A375_A2058_PLX4720_Processed.csv")


d <- subset(d, d$CellLine=="A2058")
write.csv(d, file="20150227 A2058_PLX4720_Processed.csv")
#=======================================================================================
s1 <- subset(d, d$CellLine=="A2058")
ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=conc)) +
  theme_bw()+ 
  geom_smooth(span=.15, aes(group=1),data=s1, method = "loess", size=.5, alpha=0.8, col="blue")+ 
  scale_colour_manual(values=rainbow(1)) + ylim(-2, 6)+
  theme(legend.position="none") +
  ggtitle("")+ labs(x="", y="") +
  theme(axis.text=element_text(size=12))+theme(text = element_text(size=12))

             


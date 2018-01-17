#=======================================================================================
#
# Analysis of WM Cell Lines in response to PLX4720 varying concentrations
#
#=======================================================================================
# Set the working directory, read the file, and pre-process the data
source("X:\\Users\\paudelbb\\Paper_1\\Used R Codes\\drugDose.R")
dd <- "X:/Data/CELLAVISTA/Bishal Paudel/Subclones_Bishal/LongTerm"
setwd(dd)
expName <- "20160326 SKMEL19 + 28 + PLX4720.csv"
source("X:\\Users\\paudelbb\\Paper_1\\Used R Codes\\loadCellCount.R")
d <- loadCellCount(dd, 30)

#=======================================================================================
# Assign the cell type, different concentration, drugs etc 
cell <- c('SKMEL19','SKMEL28')
d$CellLine <- rep(cell, each=15)
d$drug <- "plx"
plx <- drugDose("plx", 32, 4, 2)
d$conc <- rep(plx, each=3)
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
d$Date <- "20160326"
d <- d[, c("Date","Row","Well","CellLine","Cell.Nucleus", "Time", "drug","conc", "rep", "l2", "nl2" )]
write.csv(d, file="X:/Data/CELLAVISTA/Bishal Paudel/Subclones_Bishal/LongTerm/20160326 SKMEL19+28 + PLX4720 Processed.csv")

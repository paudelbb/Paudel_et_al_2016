#=====================================================================================================
# Figure 1: Population Level Response for different cell lines
#=====================================================================================================
#         Population level response of BRAF-mutated melanoma cell lines to PLX4720
#         CellLines used: SKMEL5, SKMEL19, SKMEL28, WM88, WM164, WM793, A375
#                 Concentrations for idling phenotype in cell lines:
#         SKMEL5(8uM), SKMEL19(8uM), SKMEL28(16uM), WM88(8uM), WM164(32uM), WM793(16uM), A375(8uM)
#=====================================================================================================
#=====================================================================================================

library(gplots)
require(ggplot2)
setwd("~/Paudel_et_al_2016/SourceData") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") 
fig = "Fig1" # Which figure is this in paper?
output = paste0(output, "/", fig)
#=====================================================================================================
# Read the experimental data
#=====================================================================================================

data <- read.csv("Population_Response_CellLines.csv", header=T, sep=",")
data = subset(data, data$Time <=350)

#=====================================================================================================
# A375 Population Level Response
#=====================================================================================================
cell <- "A375"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


dt = subset(s1, s1$Date==unique(s1$Date)[1])
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# SKMEL19 Population Level Response
#=====================================================================================================

cell <- "SKMEL19"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)

dt = subset(s1)
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# SKMEL28 Population Level Response
#=====================================================================================================

cell <- "SKMEL28"
cnc <- 16
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


dt = subset(s1)
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# WM88 Population Level Response
#=====================================================================================================

cell <- "WM88"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


dt = subset(s1)#, s1$Date==unique(s1$Date)[2])
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# WM164 Population Level Response
#=====================================================================================================

cell <- "WM164"
cnc <- 32
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)


dt = subset(s1)
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# WM793 Population Level Response
#=====================================================================================================

cell <- "WM793"
cnc <- 32
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)

dt = subset(s1)
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
# SKMEL5 Population Level Response
#=====================================================================================================

cell <- "SKMEL5"
cnc <- 8
s1 <- subset(data, data$CellLine==cell & data$conc==cnc)
s1$Date <- as.character(s1$Date)

ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)

dt = subset(s1)
y.loess <- loess(y ~ x, span=0.25, data.frame(x=dt$Time, y=dt$nl2))
y.predict <- predict(y.loess, data.frame(x=dt$Time))
plot(dt$Time,y.predict, col="red", lwd=1.)
infl <- c(FALSE, diff(diff(y.predict)>=0)!=0)
#points(dt$Time[infl ], y.predict[infl ], col="blue", lwd=3.0)
a <- c(dt$Time[infl])
unique(a)
#=====================================================================================================
#=====================================================================================================



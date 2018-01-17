#=====================================================================================================================
# Analysis of A2058 response in PLX4720
#=====================================================================================================================
# Set the working directory, read the file, and pre-process the data
dd <- "/Users/paudelbb/Paudel_et_al_2016/LongTerm"
setwd(dd)
require(gplots)
require(ggplot2)
output <- "/Users/paudelbb/Paudel_et_al_2016/LongTerm/Figs"
#=====================================================================================================================
data <- read.csv("20150227 A2058_PLX4720_Processed.csv", header=T, sep=",")
cell <- unique(data$CellLine); 
dr <- "PLX4720"
s1 <- data
cnc = "8"
#=====================================================================================================================
ggplot(data = s1, aes(x=s1$Time, y=s1$nl2, col=Date))+ 
  theme_bw()+geom_smooth(span=.25, aes(group=1), method = "loess", size=.5, alpha=0.6, col="blue")+ 
  scale_colour_manual(values="blue") + ylim(-3,3)+ xlim(0,360)+
  theme(legend.position="none") + theme(axis.text=element_text(size=12)) +
  theme(text = element_text(size=12)) + ggtitle(paste0()) + labs(x="", y="") +
  ggsave(paste0(cell, " + ", cnc, "μΜ.pdf"), path=output, width=3, height=3)
#=====================================================================================================================

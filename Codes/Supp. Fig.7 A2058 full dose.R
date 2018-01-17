#=====================================================================================================
# Load the required libraries, working directory and define output
#=====================================================================================================
library(gplots)
require(ggplot2)
require(Hmisc)
source("~//Paudel et al. 2017 Paper 1//SourceCode//summarySE.R") # Summarize function to get the summary statistics;
setwd("~/Paudel et al. 2017 Paper 1/SourceData") # Set the working directory to the folder with data;
# output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData"))), "Figures") # Set the output for plots
# fig = "Supp.Fig8" # Which figure is this in paper?
output = "/Users/paudelbb/Paudel_et_al_2016/Figures/Supp.Fig8"

#=====================================================================================================
#=====================================================================================================
# Read the experimental data
# Read the clonal fractional proliferation assay data on different melanoma celllines, 
# There are different concentrations of PLX4720 used.
#=====================================================================================================
files = dir()
popknh = load("ResponseData2016-10.Rdata")
popknh = Lizard
popknh = subset(popknh, popknh$Treatment=="PLX4720")

d = subset(popknh, popknh$Cell.Line=="A2058")

d1 = subset(d, d$Date==unique(d$Date)[3])
source("~/Paudel_et_al_2016/Codes/summarySE.R")
# d <- s1[s1$conc!=0.25, ]
d = d1
df <- summarySE(d, measurevar="nl2", groupvars=c("Time", "Conc"))
d <- df

#=====================================================================================================
# Graph function
#=====================================================================================================

par(ps = 12, cex = 1, cex.axis = 1)
graph <- function(time, nl2, conc, title, sub)
{
  y_min = -1
  y_max = 5.5
  plot(  time, nl2, 
         xlim=c(0,(max(d$Time)+35)), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub)
  with (
    data = d, expr = errbar(Time, nl2, nl2+se, nl2-se, add=T, pch=21, cap=.002))
  for(tc in unique(conc))
  {
    ccc <- conc[tc==conc]
    lines(time[tc==conc], nl2[tc==conc])
    text((max(d$Time)-0), nl2[tc==conc][length(ccc)], label=paste0(tc,'µM'), pos=4, cex=1.0)
  }
}

d$Conc = d$Conc*1000000
#=====================================================================================================
# Plot the full dose range in SKMEL5
#=====================================================================================================

pdf(paste0(output, "/A2058_full_PLX4720_Range.pdf"), height=6, width=4)
graph(d$Time, d$nl2, d$Conc, title="", sub="")
dev.off()
#=====================================================================================================

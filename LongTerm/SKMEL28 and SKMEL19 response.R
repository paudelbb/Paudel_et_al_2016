#=====================================================================================================================
# Analysis of WM lines to PLX4720 Response
#=====================================================================================================================
# Set the working directory, read the file, and pre-process the data
dd <- "/Users/paudelbb/Paudel_et_al_2016/LongTerm"
setwd(dd)
files <- grep('Processed', dir(), value=T)
files
#=====================================================================================================================
d <- read.csv("20160326 SKMEL19+28 + PLX4720 Processed.csv", header=T, sep=",")
dn <- data.frame()
for(l in unique(d$CellLine)){
  s1 <- subset(d, d$CellLine==l)
  for(i in unique(s1$conc)){
    temp <- subset(s1, s1$conc==i)
    reps <- seq(1:(length(temp$Time)/length(unique(temp$Time))))
    temp$rep <- rep(reps, each=length(unique(temp$Time)))
    dn <- rbind(dn, temp)
  }
}
d <- subset(dn, dn$Time <= 370)
unique(d$CellLine)
#=====================================================================================================================
cellchoice <- "SKMEL28"
s1 <- subset(d, d$CellLine==cellchoice)
s1 <- subset(s1, s1$Time<=200)
d <- s1

#=====================================================================================================================
# Get the fold-inhibition compared to DMSO control at different time points for IC50 plots
#=====================================================================================================================
times <- unique(d$Time)
listconc <- subset(unique(d$conc), unique(d$conc)!=0)
# X is control wells, X1 is initial value, X2 is timepoint value
x1 <- subset(d, Time == 0  & (conc == 0))
x1 <- aggregate(x1$Cell.Nucleus, by=list(Date=x1$Date, Treatment=x1$drug, CellLine=x1$CellLine, Time=x1$Time), mean)
# Get the control count for timepoint
dn <- data.frame()
for(i in 1:length(times)){
  x2 <- subset(d, Time == times[i]  & (conc == 0))
  x2 <- aggregate(x2$Cell.Nucleus, by=list(Date=x2$Date, Treatment=x2$drug, CellLine=x2$CellLine, Time = x2$Time), mean)
  x2$denom <- log2(x2$x) - log2(x1$x)
  dn <- rbind(dn, x2)
}
control <- dn
# Get the control count for timepoint
y1 <- subset(d, Time==0 & conc!=0)
dn <- data.frame()
for(i in 1:length(times)){
  for(j in unique(listconc)){
    yn <- subset(y1, conc==j)
    y2 <- subset(d, Time == times[i]  & conc == j)
    y2$numer <- log2(y2$Cell.Nucleus) - log2(mean(yn$Cell.Nucleus))
    dn <- rbind(dn, y2)
  }
}
treatment <- dn

dn <- data.frame()
for(i in unique(treatment$conc)){
  t <- subset(treatment, treatment$conc==i)
  for(j in unique(t$Time)){
    y <- subset(t, t$Time==j)
    x <- subset(control, control$Time==j)
    y$fold.inhib <- y$numer/x$denom
    y <- y[,c("Date", "Well", "CellLine", "Time", "Cell.Nucleus","drug", "conc", "rep", "fold.inhib")]
    dn <- rbind(dn, y)
  }
}

# Make numeric (This is a nasty fix to R treating concentration a factor due to the presence of "DMSO")
dn$conc <- as.numeric(as.character(dn$conc))
#=====================================================================================================================
dn$fold.inhib[dn$fold.inhib ==Inf | dn$fold.inhib ==-Inf] <- 1
# data <- subset(dn, dn$conc==16)
# data <- subset(data, data$Time<=100)
# data <- subset(data, data$Time!=0)
# plot(fold.inhib~Time, data=data)
#=====================================================================================================================
# model output is %inhibition curve (at each drug concentration)
# x = drug conc (on log scale)
# max = maximum effect of drug
# min = minimum effect of drug
# s = strength of drug effect (Hill coefficient)
#=====================================================================================================================
IC50.model	<-	function(x, IC50, max, min, s)
{
  x    <- 10^x
  IC50 <- 10^IC50
  #print(paste(IC50, max, min, s))
  (
    ( (max-min)*(x/IC50)^s ) /
      ( 1 + (x/IC50)^s )
  ) + min
}
#=====================================================================================================================

fit.plot <- function(data, title, xlab="Log10 Nanomolar", guess=1)
{
  plot(log10(data$conc), data$fold.inhib, ylim=c(-1.2,1.2), ylab="Fold Inhibition", xlab=xlab, main=title, xlim=c(-4, 4))
  abline(h=0.5, col='black', lty=2)
  
  fit <- NULL
  
  try({
    fit <-  nls(fold.inhib ~ IC50.model(x=log10(conc), IC50, max, 1.0, strength),
                data = data,
                start = list(   IC50 = guess, strength=1,
                                max=min(data$fold.inhib)),
                algorithm="port"
    )
    co <- coef(summary(fit))[,1]
    GI50 <- uniroot(function(x) IC50.model(x,co[1], co[3], 1.0, co[2])-0.5, c(-6, 6))$root
    cat("GI50=",GI50,'\n')
    attr(fit, "GI50") <- GI50
    
    x <- -200:500/100
    lines(x, IC50.model(x, co[1], co[3], co[4], co[2]))
    abline(v=co[1], col='blue')
    abline(v=GI50, col='red')
    
    curve(IC50.model(x,co[1], co[3], 1.0, co[2]), add=TRUE, col='gray', lwd=2)
  })
  
  fit
}
#=====================================================================================================================
#=====================================================================================================================
time <- unique(dn$Time)
time
data <- subset(dn, dn$Time==time[17])
fit.plot(data, title = cellchoice)
#=====================================================================================================================
#=====================================================================================================================
source("/Users/paudelbb/Paudel_et_al_2016/Codes/summarySE.R")
df <- summarySE(dn, measurevar="fold.inhib", groupvars=c("Date", "CellLine", "Time", "drug", "conc"))

#=====================================================================================================
# Graph functions to get population proliferaion dynamics with error bar
# Takes summary output that has the mean +-se or sd as desired;
#=====================================================================================================
require(Hmisc)
graph1 <- function(data, time, nl2, CellLine, title, sub, col=col, xmin, xmax, ymin, ymax)
{
  y_min =  ymin
  y_max =  ymax
  x_min =  xmin
  x_max =  xmax
  d = data
  plot(  time, nl2, 
         xlim=c(x_min,x_max), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub, col=col)
  with (
    data = d, expr = errbar(Time, nl2, nl2+sd, nl2-sd, add=T, pch=16, cap=.002, col=col))
  for(tc in unique(CellLine))
  {
    ccc <- CellLine[tc==CellLine]
    lines(time[tc==CellLine], nl2[tc==CellLine], col= col, lwd=1.0, lty=1)
    text((max(d$Time)-0), nl2[tc==CellLine][length(ccc)], label=paste0(tc), pos=4, cex=0.8)
  }
}
#=====================================================================================================
graph1(df, df$Time, df$fold.inhib, df$conc, title="", sub="", col="red", xmin=0, xmax=200, ymin=-1, ymax=2)
require(ggplot2); require(gplots)
df <- subset(df, df$conc==4)
mybarcol <- "gray20"
mp <- barplot2(df$fold.inhib, beside=TRUE,
               col=grey.colors(10), ylab="",
               main = "", font.main = 5,
               sub = "", col.sub = mybarcol, ylim=c(-1.5, 2.0),
               cex.names = 0.8, plot.ci = TRUE, ci.l = df$fold.inhib-df$se, ci.u = df$fold.inhib+df$se,
               plot.grid = FALSE, cex=0.8)
group <- c(as.character(unique(df$Time)))
group <- noquote(group)
labels <- paste0(group)
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
#=====================================================================================================
#=====================================================================================================




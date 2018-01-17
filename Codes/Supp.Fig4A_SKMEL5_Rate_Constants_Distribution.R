# ===========================================================================================================
# Figure 4A: SKMEL5 distribution of rate constants
# ===========================================================================================================
library(gplots)
require(ggplot2)
require(FME)
require(deSolve)
library(devtools)
library(ggbiplot)

setwd("~/Paudel_et_al_2016/SourceData/Model_Fitting") # Set the working directory;
# Set the output for plots
output <- paste0(substr(getwd(), 1, (nchar(getwd())-nchar("SourceData/Model_Fitting"))), "Figures") 
fig = "Supp.Fig4" # Which figure is this in paper?
output = paste0(output, "/", fig)
# ===========================================================================================================
# Get the MCMC estimated parameters
# ===========================================================================================================

set.seed(12345)
files <- grep(".RData", dir(), value=T)
chain1 <- load("MCMC_Subclones_Mix_150000_Run.RData")
dim(MCMC1$pars)
MCMC1$pars <- MCMC1$pars[-(1:dim(MCMC1$pars)[1]/2),]
rates <- data.frame(MCMC1$pars[,c(1,2,3,4)])

colnames(rates) = c("k_rs", "k_sr", "k_es", "k_se")
# ===========================================================================================================
# Plot the rate constants distributions
# ===========================================================================================================
group = names(rates)
group <- noquote(group)
labels <- group
x = c(1,2,3,4)

pdf(paste0(output, "/Boxplot_of_SKMEL5_MCMC_rate_constants.pdf"), width=3.2, height=3.2)
boxplot(rates, notch=T, col=grey.colors(4), xaxt="n")
text(x, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
dev.off()
# ===========================================================================================================
# ===========================================================================================================
boxplot(log10(1/rates), outline = F)

y1 = sin(sample(log10(1/rates$k_ab), 10000))
y2 = sin(sample(log10(1/rates$k_ba), 10000))
y3 = sin(sample(log10(1/rates$k_cb), 10000))
y4 = sin(sample(log10(1/rates$k_bc), 10000))
dt = data.frame(y1=y1, y2=y2, y3=y3, y4=y4)
dt = data.matrix(dt)
# boxplot(log10(dt), outline=F)
# boxplot(dt, outline = F, notch = T)

# ===========================================================================================================
# ===========================================================================================================

diff1 = max(rates$k_ab) - min(rates$k_ab)
diff2 = max(rates$k_ba) - min(rates$k_ba)
diff3 = max(rates$k_cb) - min(rates$k_cb)
diff4 = max(rates$k_bc) - min(rates$k_bc)

# ===========================================================================================================
library(preprocessCore)
nd = normalize.quantiles(dt,copy=TRUE)
nd = data.frame(nd)
boxplot(nd, outline = F, notch = T)
# ===========================================================================================================
# plot(y1, type="l", xlim=c(0,1000))
# lines(y2, col="red")
# lines(y3, col="blue")
# lines(y4, col="green")

# ===========================================================================================================
# ===========================================================================================================

k1 = sample(log10(1/rates$k_ab), 10)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k1)){
  dc = 1/2*k1[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k1[i])
  dn = rbind(dn, dc)
}
df1 = dn
plot(y~x, data=df1, ylim=c(-100, 150))
#
k2 = sample(log10(1/rates$k_ba), 10)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k2)){
  dc = 1/2*k2[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k2[i])
  dn = rbind(dn, dc)
}
df2 = dn
df2$x = df2$x+20
plot(y~x, data=df2, ylim=c(-100, 150))
# ===========================================================================================================
# ===========================================================================================================
k3 = sample(log10(1/rates$k_bc), 10)
x1 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k3)){
  dc = 1/2*k3[i]*x1^2
  dc = data.frame(x = x1, y=dc, k=k3[i])
  dn = rbind(dn, dc)
}
df3 = dn
df3$x = df3$x + 20
plot(y~x, data=df3, ylim=c(-100, 250))
#
k4 = sample(log10(1/rates$k_cb), 10)
x2 = seq(-10, 10)
dn = data.frame()
for(i in 1:length(k4)){
  dc = 1/2*k4[i]*x2^2
  dc = data.frame(x = x2, y=dc, k=k4[i])
  dn = rbind(dn, dc)
}
df4 = dn
df4$x = df4$x +30 + 10
plot(y~x, data=df4, ylim=c(-100, 250))

# ===========================================================================================================
# Quantile normalize the Well B
# ===========================================================================================================
dty = data.frame(b1=df2$y, b2=df3$y)
library(preprocessCore)
dt = data.matrix(dty)
nd = normalize.quantiles(dt,copy=TRUE)
nd = data.frame(nd)


df1
df2$y = nd$X1
df3$y = nd$X2
df4
# ===========================================================================================================
dn = data.frame()
for(i in 1:length(unique(df1$k))){
  temp1 = subset(df1, df1$k==unique(df1$k)[i])
  for(j in 1:length(unique(df2$k))){
    temp2 = subset(df2, df2$k==unique(df2$k)[j])
    diff = max(temp2$y) - max(temp1$y)
    temp2$y = temp2$y-diff
    dc = rbind(temp1,temp2)
    dn = rbind(dn, dc)
  }
}
uid = seq(1:(length(dn$x)/42))
dn$id = rep(uid, each=42)

dt1 = dn
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-100, 250))
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}


# ===========================================================================================================
# ===========================================================================================================
dn = data.frame()
for(i in 1:length(unique(df3$k))){
  temp1 = subset(df3, df3$k==unique(df3$k)[i])
  for(j in 1:length(unique(df4$k))){
    temp2 = subset(df4, df4$k==unique(df4$k)[j])
    diff = max(temp2$y) - max(temp1$y)
    temp2$y = temp2$y - diff
    dc = rbind(temp1,temp2)
    dn = rbind(dn, dc)
  }
}
uid = seq(1:(length(dn$x)/42))
dn$id = rep(uid, each=42)

dt2 = dn
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-100, 250))
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}


# ===========================================================================================================
# ===========================================================================================================


dn = data.frame()
for(i in 1:length(unique(dt1$id))){
  temp1 = subset(dt1, dt1$id==unique(dt1$id)[i])
  for(j in 1:length(unique(dt2$id))){
    temp2 = subset(dt2, dt2$id==unique(dt2$id)[j])
    diff = max(temp2$y) - max(temp1$y)
    temp2$y = temp2$y - diff
    dc = rbind(temp1,temp2)
    dn = rbind(dn, dc)
  }
}
uid = seq(1:(length(dn$x)/42))
dn$id = rep(uid, each=42)
df = dn
cl1 <- rainbow(length(unique(df$id)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(y~x, data=df, type="n", ylab="", xlab="", ylim=c(-250, 250), xaxt="n")
for(w in unique(df$id)){
  temp <- subset(df, df$id==w)
  lines(temp$x, temp$y, col=cl1[w])
}


#=====================================================================================================
# Graph functions to get population proliferaion dynamics with error bar
# Takes summary output that has the mean +-se or sd as desired;
#=====================================================================================================
require(Hmisc)
graph2 <- function(data, time, nl2, CellLine, title, sub, col=col, xmin, xmax, ymin, ymax)
{
  y_min =  ymin
  y_max =  ymax
  x_min =  xmin
  x_max =  xmax
  d = data
  points(  time, nl2, 
         xlim=c(xmin,xmax), ylim=c(y_min, y_max),
         xlab='',ylab='',
         main=title, sub=sub, col=col)
  with (
    data = d, expr = errbar(Time, nl2, nl2+sd, nl2-sd, add=T, pch=16, cap=.002, col=col))
  for(tc in unique(CellLine))
  {
    ccc <- CellLine[tc==CellLine]
    lines(time[tc==CellLine], nl2[tc==CellLine], col= col, lwd=1.0, lty=2)
    text((max(d$Time)-0), nl2[tc==CellLine][length(ccc)], label=paste0("SC", substr(tc,16, 17),''), pos=4, cex=0.8)
  }
}
#=====================================================================================================
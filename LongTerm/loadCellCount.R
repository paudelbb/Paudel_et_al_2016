loadCellCount <- function(folderWithData=getwd(), wellCount)
{
startDir <- getwd()
if(folderWithData!=startDir) 	setwd(folderWithData)
d <- read.csv(expName, header=T, sep=",")
d <- d[, c("Row", "Column", "Cell.Nucleus")]
d$Row <- as.character(d$Row)
d$Column <- as.character(d$Column)
idx <- grep('Timespan', d[,2])
times <- as.numeric(d[,1][idx])
time.pt <- 1:length(times)
for(tp in time.pt)
	{	
        temp			<-	d[(idx[tp]+1):(idx[tp]+wellCount),]
		temp$Time		<-	times[tp]
		ifelse(tp==time.pt[1], a <- temp, a <- rbind(a,temp))
	}
rownames(a)	<-	seq(nrow(a))
a[nchar(a$Column)==1,'Column']	<-	paste0('0',a[nchar(a$Column)==1,'Column'])
a$Well		<-	paste0(a$Row,a$Column)
a$Column	<-	as.integer(a$Column)
a$l2		<-	log2(a$Cell.Nucleus)
a
}
args<-commandArgs(TRUE)
input=args[1]
output=args[2]

Cscore<-read.delim(input,header=F)

cutoff<-NULL

# 99.9%~0.1%

for (i in 1:999)
{
	cutoff<-rbind(cutoff,sort(t(Cscore),decreasing=TRUE)[ceiling(nrow(Cscore)*(1-i/1000))])
}


write(t(cutoff),file=output, ncolumns=1)

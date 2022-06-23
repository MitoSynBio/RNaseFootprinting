#!/usr/bin/Rscript

args<-commandArgs(TRUE)
input=args[1]
output=args[2]

Cscore <- read.delim(input,header=F)

# Shuffle last column (Cscore)

shuf <- data.frame(A=Cscore$V1, B=Cscore$V2, C=Cscore$V3, D=sample(Cscore$V4))

write.table(shuf, output, col.names=F, row.names=F, sep="\t", quote=F)



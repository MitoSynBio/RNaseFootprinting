#!/usr/bin/Rscript

args<-commandArgs(TRUE)
input1=args[1]
input2=args[2]
output=args[3]

exp <- read.delim(input1, header=T)
null <- read.delim(input2, header=T)


## x_plus_real
count <- c()
for (i in 1:length(exp[,10])) {
  count[i] <- sum(exp[,10] >= exp[i,10])
}
exp$x_plus_real <- c(count)

## x_plus_null
count2 <- c()
for (i in 1:length(exp[,10])) {
  count2[i] <- sum(null[,8] >= exp[i,10])
}
exp$x_plus_null <- c(count2)

## FDR
FDR_count <- c()
for (i in 1:length(exp[,10])) {
  FDR_count[i] <- exp[i,12]/exp[i,11]
}
exp$FDR <- c(FDR_count)

## FDR_norm
norm_count <- c()
real_len <- length(exp[,10])
null_len <- length(null[,8])

for (i in 1:length(exp[,10])) {
  norm_count[i] <- exp[i,13]*(real_len/null_len)
}
exp$FDR_norm <- c(norm_count)

write.table(exp, output, col.names=T, row.names=F, sep="\t", quote=F)



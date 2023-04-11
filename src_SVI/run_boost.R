library(utils)
library(stats)
library(grDevices)


args = commandArgs()

if (length(args)==0) {
  stop("not enough input", call.=FALSE)
}

count_f <- args[4]
meta_f <- args[5]
species <- args[6]
out_f <- args[7]

colData <- read.csv(meta_f, stringsAsFactors=FALSE, row.names=1, check.names=F)
counts <- read.csv(count_f, row.names=1, check.names=F, stringsAsFactors=FALSE)

res <- boost.gp(Y = counts, loc = colData)
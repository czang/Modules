#!/bin/sh

INFILE=$1_peak_lengths.bed
OUTFILE=$1_peak_lengths_hist.pdf

/usr/local/bin/R --no-save -q << EOF

data <- read.table("$INFILE", header=T)
Llist <- data[,4]

pdf("$OUTFILE")
hist(Llist)
dev.off()

print(c("Median","$1", median(Llist)))
EOF

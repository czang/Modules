#!/bin/sh

MOTIF=$1
INFILE=$MOTIF
OUTFILE=utest_$MOTIF
ALLFILE=../motif_hit_count.txt

COLUM=3

ALLCOUNT=646923
COUNT=5000

/usr/bin/R64 --no-save -q << EOF

b <- read.table("$ALLFILE")[,2]
x <- read.table("$INFILE")[,$COLUM]
l <- length(x)
down.p <- wilcox.test(x,b,alternative="less")\$p.value
up.p <-  wilcox.test(x,b,alternative="greater")\$p.value
r <- cbind("$MOTIF",l,down.p,up.p)
write.table(r,"$OUTFILE",sep='\t',quote=F,row.names=F,col.names=F)

EOF

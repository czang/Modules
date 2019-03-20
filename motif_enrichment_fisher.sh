#!/bin/sh

INFILE=$1
OUTFILE=$2
ALLFILE=../motif_hit_count.txt
ALLCOUNT=646923
COUNT=5000

/usr/bin/R64 --no-save -q << EOF

g <- read.table("$ALLFILE")[,2]
rawgata <- read.table("$INFILE")[,2:5]
rownames(rawgata) <- read.table("$INFILE")[,1]
ffp <- function(x){fisher.test(matrix(c(x,$COUNT,$ALLCOUNT),nrow=2),alternative="greater")\$p.value}
p.1 <- apply(cbind(rawgata[,1],g),1,ffp)
p.2 <- apply(cbind(rawgata[,2],g),1,ffp)
p.3 <- apply(cbind(rawgata[,3],g),1,ffp)
p.4 <- apply(cbind(rawgata[,4],g),1,ffp)
r <- rawgata/$COUNT
gr <- g/$ALLCOUNT
f <- log2(r/gr)
final <- as.data.frame(cbind(f[,1],p.1,f[,2],p.2,f[,3],p.3,f[,4],p.4))
rownames(final) <- rownames(rawgata)
write.table(final,"$OUTFILE",sep='\t',quote=F)

EOF

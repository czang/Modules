#!/bin/sh

if [ $# -lt 2 ]; then
        echo "Need command line arguments for file names"
        exit 1
fi


INFILE=$1
OUTFILE=$2

#lanes=(1 2 3 4 5 6 7)
lanes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
DEST_DIR=.

for ((i=0;i<${#lanes[@]};i+=1)); do
	CHROM=chr${lanes[$i]}
	grep "$CHROM[[:space:]]" $INFILE > $CHROM.tmp
	sort -g -k 2 $CHROM.tmp > $CHROM.sort
	#cat $CHROM.sort >> $OUTFILE
	#rm $CHROM.*
done

cat chr*.sort > $OUTFILE
rm chr*.sort
rm chr*.tmp
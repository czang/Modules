#!/bin/sh
if [ $# -lt 1 ]; then
        echo "Need a command line argument for library name, e.g., H3K4me3"
        exit 1
fi

INPUT=$1.bed
DIR=.
RESOLUTION=10
CHROM=chr1

grep ${CHROM}[[:space:]] $DIR/$INPUT > temp1.bed
echo "Total tag count in $CHROM is"
wc -l temp1.bed
grep + temp1.bed > temp_plus.bed
grep - temp2.bed > temp_minus.bed

echo "Positive"
python /home/zang/Modules/make-graph-file.py -f temp_plus.bed -c $CHROM -w $RESOLUTION -i 0 -o temp_plus.graph

echo "Negative"
python /home/zang/Modules/make-graph-file.py -f temp_minus.bed -c $CHROM -w $RESOLUTION -i 0 -o temp_minus.graph

echo "Calculating correlation"
python /home/zang/Modules/calculate_cross_correlation_long_range.py -s hg18 -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $1-$CHROM-tag-correlation.txt

rm temp*
sh make-correlation-plot.sh $1-$CHROM-tag-correlation
ps2pdf $1-$CHROM-tag-correlation.ps

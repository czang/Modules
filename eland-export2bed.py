#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


## Important: Right now this script will only pick up reads that align to
## traditional chromosome files.  That means this will miss any reads that
## match to files such as splice junction files.


## TO DO: Standardize a method to treat reads aligning to splice junctions.


"""
According to the Illumina documentation, the most usable output file from
using eland_extended is:

s_N_sorted.txt


This script is to make our usual bed files from these files
"""


"""
example lines from s_N_sorted.txt

HWI-EAS434      246     2       40      992     335     0       1       GCTAGGAGTTGACTACAT      abb`Tababbb_bbbbbb      chr10.fa          3445460  F       18      11

HWI-EAS434      246     2       17      384     916     0       1       GTCCACGATCCTGTAGGC      aaabbbaa`bbaaYa[D^      chr14.fa          97786820 F       16T1    4


machine - run - lane - tile - x-coord - y-coord - index - read# - read -
quality string - match chrom - match position (1 start) - strand - match
descriptor - single-read alignment score


~~~~~~~~~~~~~ match descriptor ~~~~~~~~~~~~~~~

-- for 18 bp read, 18 indicates perfect 18 bp match

-- 16T1 = indicates a substitution of T at the 16th position


~~~~~~~~ single-read alignment score ~~~~~~~~~

Alignment score of a single-read mtch, or for a paired read, alignment score
of a read if it were treated as a single read.  Blank if no match found; any
score less than 4 should be considered as aligned to a repeat.


(see pg 119 of manual for full description)

"""


def getBeds(sfile):
    """ BED6 = chrom + start + end + name + score + strand """
    beds = [];
    infile = open(sfile);
    for line in infile:
        line = line.strip();
        sline = line.split();
        chrom = re.sub(".fa", "", sline[10]);
        
        if re.match ("c", chrom):
            """ putting this in to check to make sure the read is aligning to
            a proper chromosome file -- it could also be aligning to a splice
            junction file, and then has to be treated separately"""
            
            chrom = re.sub("c", "chr", chrom)
            start = atoi(sline[11]) - 1; ## -1 to deal with Eland starting at 1 and UCSC starting at 0
            name = sline[13];  ## name is match descriptor
            score = 0;
            if sline[14]:
                score = atoi(sline[14]);
            strand = '';
            if re.match("F", sline[12]):
                strand = "+";
            elif re.match("R", sline[12]):
                strand = "-";
            end = start + len(sline[8]); ## should really be length of alignment
            bed = chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\t" + str(score) + "\t" + strand;
            beds.append(bed);
    return beds;



def printOutBeds(beds, outfilename):
    outfile = open(outfilename, 'w');
    for b in beds:
        outline = b + "\n";
        outfile.write(outline);
    outfile.close();


def main(argv):
    parser = OptionParser()
    parser.add_option("-s", "--eland_sorted_out", action="store", type="string",
                      dest="sfile", help="eland sorted file",
                      metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)

    beds = getBeds(opt.sfile);
    printOutBeds(beds, opt.outfile);


    

if __name__ == "__main__":
    main(sys.argv)


        

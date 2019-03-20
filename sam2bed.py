#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser


"""
Convert the SAM format output of tophat to BED.
"""

# IMPORTANT: this will only work for single-end reads

def getStrand(flag):
    """
    Convert flag to strand
    """
    strand = "+"
    if (flag & (0x0010)):
        strand = "-"		
    return strand
    
def doConversion(File, outfilename):
    infile = open(File);
    outfile = open(outfilename, 'w');
    for line in infile:
        line = line.strip();
        sline = line.split();
        chrom = sline[2];
        start = atoi(sline[3]) - 1 ## make 0-based
        end = start + len(sline[9]);
        name = sline[0];
        strand = getStrand(atoi(sline[1]));
        score = 0; ## place holder
        outline = chrom + "\t" + str(start) + "\t" + str(end) + "\t" + \
                  name + "\t" + str(score) + "\t" + strand + "\n";
        outfile.write(outline);
    outfile.close();

def main(argv):
    parser = OptionParser()
    parser.add_option("-s", "--samfile", action="store", type="string",
                      dest="samfile", help="sam file", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="outfile name", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)
    doConversion(opt.samfile, opt.outfile);


if __name__ == "__main__":
    main(sys.argv)


        

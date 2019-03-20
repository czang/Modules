#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# schonesde@mail.nih.gov).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import BED;


"""
A script to get sequence for bed files.

Take in a bed file and get a directory of nib files from the command
line.  Get the sequences corresponding to the bed coordinates using
nibFrag and print out the sequences in fasta format.
"""

def getSequences(bed_vals, nibdir):

    unique_seqs = {};
    for chrom in bed_vals.keys():
        coords = bed_vals[chrom];
        for c in coords:
            strand = "+";
            nibFile = chrom + ".nib";
            nibFile = os.path.join(nibdir,nibFile);
            tmpFile = '/tmp/temp%s.%s' % (os.getpid(), c.start)
            try:
                if os.system('/home/dschones/bin/x86_64/nibFrag %s %s %s %s %s' %
                             (nibFile, c.start, c.end, strand, tmpFile)): raise
            except: sys.stderr.write("nibFile failed\n")
            infile = open(tmpFile, 'r');
            header = infile.readline();
            seq = '';
            for line in infile:
                line = line.strip();
                seq = seq + upper(line);
            seq_id = ">" + chrom + ":" + str(c.start) + "-" + str(c.end);
            unique_seqs[seq_id] = seq;
            
            ## clean up the junk
            try:
                if os.remove('%s' %
                             (tmpFile)): raise
            except: sys.stderr.write("clean up failed\n")

            """ then add reverse complement """
            '''
	    tmpFile = '/tmp/temp%s.%s.rc' % (os.getpid(), c.start)
            fake_minus = "-";
            try:
                if os.system('/home/dschones/bin/x86_64/nibFrag %s %s %s %s %s' %
                             (nibFile, c.start, c.end, fake_minus, tmpFile)): raise
            except: sys.stderr.write("nibFile failed\n")
            infile = open(tmpFile, 'r');
            header = infile.readline();
            seq = '';
            for line in infile:
                line = line.strip();
                seq = seq + upper(line);
            seq_id = ">" + chrom + ":" + str(c.start) + "-" + str(c.end) + "." + "rc";
            unique_seqs[seq_id] = seq;
            '''
            ##clean up junk
            try:
                if os.remove('%s' %
                             (tmpFile)): raise
            except: sys.stderr.write("clean up failed\n")
    return unique_seqs;


def printOut(seqs, outfilename):
    outfile = open(outfilename, 'w');
    for k in seqs.keys():
        outline = k + "\n";
        outfile.write(outline);
        outline = seqs[k] + "\n";
        outfile.write(outline);
    outfile.close();


def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--bed_file", action="store", type="string",
                      dest="bedfile", help="bed file", metavar="<file>")
    parser.add_option("-d", "--nibdir", action="store", type="string",
                      dest="nibdir", help="directory containing nib files",
                      metavar="<dir>")
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species",
                      metavar="<str>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="file to print output to", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)


    bed_vals = BED.BED(opt.species, opt.bedfile, "BED3");
    seqs = getSequences(bed_vals, opt.nibdir);
    printOut(seqs, opt.outfile);


if __name__ == "__main__":
    main(sys.argv)


        

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import BED
import GenomeData
import FragSizeEstimation


"""
(1) Break bed file up into chromosomes
(2) Get the average fragment size
(3) Make summary files for loading into genome browser
"""


Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";


make_graph_file = os.path.join(Dir, "make-nr-graph-file.py");


def separateByChrom(chroms, file):
    for chrom in chroms:
        match = chrom + "[[:space:]]";
        tmpFile = chrom + ".bed";
        try:
            if os.system('%s %s %s > %s' %
                         (grep, match, file, tmpFile)): raise
        except: sys.stderr.write("grep failed at %s\n % chrom");


def getAveFragSize(chroms, species, small, large):
    ave_frag_size = 0;
    for chrom in chroms:
        bed_file = chrom + ".bed";
        frag = FragSizeEstimation.FragSize(bed_file, species);
        ave_frag_size += frag.getMode(small, large);

    ave_frag_size = int(float(ave_frag_size) / float(len(chroms)));
    return ave_frag_size;


def makeGraphFile(chroms, fragment_size):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        graph_file = chrom + ".graph";
        try:
            if os.system('%s -f %s -c %s -i %s -o %s'%
                         (make_graph_file, bed_file, chrom,
                          fragment_size, graph_file)): raise
        except: sys.stderr.write("Summary file construction failed at %s\n" % chrom)



def combineAllGraphFiles(chroms, final_out):
    """
    Combine the seperately processed chromosomes
    """
    outfile = open(final_out,'w');
    outfile.close();
    
    for chrom in chroms:
        graph_file = chrom + ".graph";
        try:
            if os.system('%s %s >> %s' %
                         (cat, graph_file, final_out)): raise
        except: sys.stderr.write("cat failed at %s\n" % chrom)
    
		

def cleanup(chroms):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        graph_file = chrom + ".graph";
        try:
            os.remove('%s' % bed_file)
            os.remove('%s' % graph_file)
        except: sys.stderr.write("clean up failed at %s\n" % chrom);



    
def main(argv):
    parser = OptionParser()
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="mm8,hg18,dm3,etc", metavar="<str>")
    parser.add_option("-b", "--bed_file", action="store", type="string",
                      dest="bedfile", help="bed file to make graph file of",
                      metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output bed summary file name",
                      metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 6:
        parser.print_help()
        sys.exit(1)

    small_size = 20;
    large_size = 500;

    if opt.species in GenomeData.species_chroms.keys():
        chroms = GenomeData.species_chroms[opt.species];
        separateByChrom(chroms, opt.bedfile);
        
        ave_frag_size = getAveFragSize(chroms, opt.species, small_size, large_size);
        makeGraphFile(chroms, frag_size);
        combineAllGraphFiles(chroms, opt.outfile);
        cleanup(chroms);
    else:
        print opt.species + " is not in the species list";
    
    

if __name__ == "__main__":
	main(sys.argv)


        

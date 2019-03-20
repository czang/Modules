#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


def get_genome_length(filename):
	infile = open(filename,'r');
	total_seq_length = 0.0;
	for line in infile:
		line = line.strip();
		sline = line.split();
		total_seq_length += atof(sline[1]);
	infile.close();
	return total_seq_length;
	
def main(argv):
    	parser = OptionParser()
    	parser.add_option("-f", "--Chrom_length_file", action="store", type="string",
                      dest="filename", help="Chrom_length_file",
                      metavar="<file>")
    
    	(opt, args) = parser.parse_args(argv)
    	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	print get_genome_length(opt.filename);


if __name__ == "__main__":
    main(sys.argv)

		

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


def get_average_tag_density (filename):
	print filename;
	infile = open(filename,'r');
	total_tag_count = 0.0;
	total_seq_length = 0.0;
	for line in infile:
		line = line.strip();
		sline = line.split();
		total_tag_count += atof(sline[1]);
		total_seq_length += atof(sline[2]);
	infile.close();
	return total_tag_count/total_seq_length;
	
def main(argv):
    	parser = OptionParser()
    	parser.add_option("-f", "--tag_count_statistics_file", action="store", type="string",
                      dest="filename", help="chrom tag_count chrom_length tag_density_in_chrom",
                      metavar="<file>")
    
    	(opt, args) = parser.parse_args(argv)
    	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	print opt.filename;	
	print get_average_tag_density(opt.filename);


if __name__ == "__main__":
    main(sys.argv)

		

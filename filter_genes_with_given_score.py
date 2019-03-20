#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import GenomeData;
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--genesummary", action="store", type="string", dest="islandsummary", metavar="<file>", help="gene summary file")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", metavar="<int>", help="column number, start from 0")
	parser.add_option("-t", "--score", action="store", type="int", dest="threshold", metavar="<int>", help="threshold")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="genes under threshold")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.islandsummary,'r');
	outfile = open(opt.out_file, 'w');
	total =0;
	for line in inputfile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > opt.column:
				flag = sline[opt.column]
				if not re.match('None', flag):
					if atoi(flag) == opt.threshold:
						total += 1;
						outfile.write('\t'.join(sline)+'\n');
	inputfile.close()
	outfile.close()
	
	print "Total number of genes with score", opt.threshold, "is",  total;
	
if __name__ == "__main__":
	main(sys.argv)

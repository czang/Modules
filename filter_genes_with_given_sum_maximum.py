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
	parser.add_option("-t", "--score", action="store", type="float", dest="threshold", metavar="<float>", help="threshold")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="genes under threshold")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.islandsummary,'r')
	outfile = open(opt.out_file, 'w')
	total = 0
	for line in inputfile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List = []
			for i in range(1,len(sline)):
				List.append(atof(sline[i]))
			if sum(List) >= opt.threshold: 
				total += 1;
				outfile.write('\t'.join(sline)+'\n');
	inputfile.close()
	outfile.close()
	
	print "Total number of genes with sum min score", opt.threshold, "is",  total;
	
if __name__ == "__main__":
	main(sys.argv)

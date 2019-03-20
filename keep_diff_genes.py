#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import separate_by_chrom


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to separate BED file into 2, according to the size of each region. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="infile", metavar="<file>", help="sample name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.infile, 'r')
	u = open(opt.out_file, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			datalist = []
			for i in range(1,len(sline)):
				datalist.append(atof(sline[i]))
			if max(datalist) >= 6.0 and max(datalist) - min(datalist) >= 1.0:
				u.write('\t'.join(sline) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
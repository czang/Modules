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
		if not re.search("s0", line):
			line = line.strip()
			sline = line.split()
			name = sline[1].strip('"')
			i = 2
			s = 1000
			while i < len(sline):
				if abs(atof(sline[i].strip('"'))) >= 0.00001:
					s = i - 1
					break
				else:
					i += 1
			u.write(name + '\t' + str(s) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
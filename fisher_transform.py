#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import stats


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to calculate fisher transformation.
z = 1/2 * ln(1 + r) / ln (1 - r)

'''

def fisher_transform(r):
	return log(1.0 + r) / log (1.0 - r) / 2.0

def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="bedfile", metavar="<file>", help="sample name")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", metavar="<int>", help="column label, start from 0", default=1)
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.bedfile, 'r')
	u = open(opt.out_file, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			r = atof(sline[opt.column])
			sline[opt.column] = str(fisher_transform(r))
			u.write('\t'.join(sline) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
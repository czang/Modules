#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to separate BED file into 2, according to the size of each region. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="infile", metavar="<file>", help="sample name")
	parser.add_option("-e", "--halfsize", action="store", type="int", dest="shift", metavar="<int>", help="half of region size")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.infile, 'r')
	u = open(opt.out_file, 'w')
	for line in f:
		if not re.search("#", line):
			line = line.strip()
			sline = line.split()
			start = atoi(sline[1])
			end = atoi(sline[2])
			center = (start+end) / 2
			sline[1] = str(center - opt.shift)
			sline[2] = str(center + opt.shift)
			u.write('\t'.join(sline) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
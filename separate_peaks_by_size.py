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
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="bed file")
	parser.add_option("-c", "--cutoff", action="store", type="int", dest="cutoff", help="cutoff", metavar="<int>")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.bedfile, 'r')
	u = open(opt.out_file+'_above_'+str(opt.cutoff), 'w')
	d = open(opt.out_file+'_below_'+str(opt.cutoff), 'w')
	for line in f:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			if int(sline[2]) - int(sline[1]) > opt.cutoff:
				u.write('\t'.join(sline)+'\n')
			else:
				d.write('\t'.join(sline)+'\n')
	f.close()
	u.close()
	d.close()


if __name__ == "__main__":
	main(sys.argv)
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
	parser.add_option("-i", "--sample", action="store", type="string", dest="infile", metavar="<file>", help="RSEG output file")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="domain output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	p = open(opt.infile, 'r')
	o = open(opt.out_file, 'w')
	for line in p:
		line = line.strip()
		sline = line.split()
		if sline[3] == 'SAMPLE-I-ENRICHED':
			o.write('\t'.join(sline)+'\n')
	p.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
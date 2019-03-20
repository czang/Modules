#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


def add_line_numbers(infile, threshold, outfile):
	f = open(infile,'r')
	o = open(outfile,'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if atoi(sline[3]) > threshold:
				o.write('\t'.join(sline)+'\n')
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-t", "--threshold", action="store", type="int", dest="threshold", help="window threshold", metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	add_line_numbers(opt.infile, opt.threshold, opt.outfile)


if __name__ == "__main__":
	main(sys.argv)
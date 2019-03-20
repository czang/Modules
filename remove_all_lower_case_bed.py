#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


def remove_all_lower_case_bed(infile, outfile):
	f = open(infile,'r')
	o = open(outfile,'w')
	for line in f:
		line = line.strip()
		if (('A' in line) or ('C' in line) or ('G' in line) or ('T' in line)):
			o.write(line + '\n')
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	remove_all_lower_case_bed(opt.infile, opt.outfile)


if __name__ == "__main__":
	main(sys.argv)
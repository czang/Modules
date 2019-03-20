#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator




def sum_column(input_file, number):
	"""
	number is the 1-based column number
	"""
	
	infile = open(input_file,'r')
	total = 0.0
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= number:
				total += atof(sline[number-1])
	infile.close()
	return total


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	print opt.input_file, sum_column(opt.input_file, opt.colum)


if __name__ == "__main__":
	main(sys.argv)

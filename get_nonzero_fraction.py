#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

#from gene_set_manipulation import *
#from associate_binary_modification_with_expression import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def get_float_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[c]))
	file.close()
	return List


def get_nonzero_fraction(gene_file, c):
	file = open(gene_file,'r')
	total = 0
	number = 0
	for line in file:
		if not re.match("#", line):
			total += 1
			line = line.strip()
			sline = line.split()
			if atof(sline[c]) > 0.01:
				number += 1
	file.close()
	fraction = float(number)/float(total)
	return fraction


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	print get_nonzero_fraction(opt.inputfile, 1)


if __name__ == "__main__":
	main(sys.argv)
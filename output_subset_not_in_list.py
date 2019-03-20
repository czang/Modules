#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="input bed file with complete data")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="subsetfile", metavar="<file>", help="bed file with ID in the 1st column")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_gene_list(opt.subsetfile,0)
	output_noinsubset_in_file (opt.inputfile, 0, List, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)

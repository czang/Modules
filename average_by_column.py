#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	datalist = get_float_list(opt.input_file, opt.colum-1)
	print average(datalist), standard_deviation(datalist)


if __name__ == "__main__":
	main(sys.argv)

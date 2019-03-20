#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--completedata", action="store", type="string", dest="infile", metavar="<file>", help="input whole data set file name")
	parser.add_option("-n", "--number", action="store", type="int", dest="number", metavar="<int>", help="number of modifications")
	parser.add_option("-p", "--pattern", action="store", type="string", dest="pattern", metavar="<str>", help="combinatorial patternin string")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	datadic = read_complete_binary_file(opt.infile)
	pattern_list = get_pattern_list_from_string(opt.pattern, opt.number)
	genelist = get_pattern_gene_set(datadic, pattern_list)
	print len(genelist)


if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import associate_binary_modification_with_expression

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input_file", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	datadic = associate_binary_modification_with_expression.read_complete_binary_file(opt.infile)
	List = associate_binary_modification_with_expression.get_modification_gene_fraction_list(datadic)
	
	f = open(opt.out_file, 'w')
	for i in List:
		f.write(str(i) + '\t')
	f.write('\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
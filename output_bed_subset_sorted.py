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
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="input bed file with complete data, 1st column must be sorted")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="subsetfile", metavar="<file>", help="bed file with ID in the 4th column, must be sorted")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_gene_list(opt.subsetfile,0)
	f = open(opt.inputfile, 'r')
	o = open(opt.out_file, 'w')
	i = 0
	for line in f:
		if not re.match('#', line):
			line = line.strip()
			sline = line.split()
			if i < len(List) and sline[3] == List[i]:
				o.write(line + '\n')
				i += 1
	f.close()
	o.close()
	

if __name__ == "__main__":
	main(sys.argv)

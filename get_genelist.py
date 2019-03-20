#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import gene_set_manipulation


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--genefile", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	List = gene_set_manipulation.get_gene_list(opt.infile, 0)
	o = open(opt.out_file, 'w')
	o.write(' '.join(List)+'\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)
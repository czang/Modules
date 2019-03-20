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
	parser.add_option("-a", "--patterndata1", action="store", type="string", dest="infile1", metavar="<file>", help="whole data set file of initial")
	parser.add_option("-b", "--patterndata2", action="store", type="string", dest="infile2", metavar="<file>", help="whole data set file of final")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	datadic1 = read_complete_binary_file(opt.infile1)
	datadic2 = read_complete_binary_file(opt.infile2)
	f = open(opt.output, 'w')
	for gene in datadic1.keys():
		if gene in datadic2.keys():
			f.write(gene)
			for i in range(0, len(datadic1[gene])):
				f.write('\t'+str(datadic2[gene][i]-datadic1[gene][i]))
			f.write('\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
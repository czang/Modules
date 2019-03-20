#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
from stats import *

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile", action="store", type="string", dest="IDfile", metavar="<file>", help="input file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	sample = get_float_list("candidate.txt", 2)
	population_list = get_float_list(opt.IDfile, 2)
	mu = mean(population_list)
	sigma = stdev(population_list)
	x = mean(sample)
	z = (x - mu)/sigma
	print opt.IDfile, z


if __name__ == "__main__":
	main(sys.argv)
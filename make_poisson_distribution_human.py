#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from species_chroms import *
from associate_island_with_genes import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	genomelength = 0.0
	for item in species_chrom_lengths['hg18']:
		genomelength += item
	average = 1000000.0/genomelength*200
	
	f = open(opt.out_file, 'w')
	for i in range(0,3000):
		f.write(str(i)+'\t'+str(poisson(i, average)*100)+'\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
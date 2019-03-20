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
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-l", "--infile2", action="store", type="string", dest="subsetfile", metavar="<file>", help="input file with genelist data")
	parser.add_option("-n", "--groupsize", action="store", type="int", dest="groupsize", metavar="<int>", help="number of genes in one group")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_gene_list(opt.subsetfile,0)
	groupnumber = len(List)/opt.groupsize
	g = 1
	i = 0
	while g < groupnumber:
		sublist = []
		for j in range(0, opt.groupsize):
			sublist.append(List[i])
			i += 1
		output_subset_in_file(opt.inputfile, 3, sublist, opt.out_file+'_'+str(g))
		g += 1
	sublist = []
	while i < len(List):
		sublist.append(List[i])
		i +=1
	output_subset_in_file(opt.inputfile, 3, sublist, opt.out_file+'_'+str(g))


if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *

def cdf(List):
	List.sort()
	cdflist = []
	cdflist.append(0.0)
	length = len(List)
	for i in range(1,100):
		j = length / 100 * i
		cdflist.append(List[j])
	cdflist.append(List[-1])
	return cdflist


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--wiggle", action="store", type="string", dest="wiggle_file", metavar="<file>", help="input wiggle file")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output wiggle file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
		
	inputlist = get_float_list(opt.wiggle_file,0)
	cdflist = cdf(inputlist)
	outfile = open(opt.out_file, 'w')
	#
	for i in range(0,101):
		outfile.write(str(cdflist[i]) + '\t' + str(i/100.0) + '\n')
	outfile.close()

	
if __name__ == "__main__":
	main(sys.argv)

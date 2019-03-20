#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *

## get BED module
#import BED;
#import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def gene_modification(datafile):
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result[sline[0]] = atof(sline[1])
	file.close()
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfiles", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-l", "--patternlist", action="store", type="string", dest="patterns", metavar="<file>", help="list of all pattern variables")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	patternlist = get_gene_list(opt.patterns, 0)
	
	complete_data = {}
	for mod in patternlist:
		complete_data[mod] = gene_modification(mod+opt.infile)
	
	f = open(opt.output,'w')
#	f.write('#gene')
#	for mod in patternlist:
#		f.write('\t' + mod)
	for gene in complete_data[patternlist[0]].keys():
#		f.write('\n'+gene)
		for mod in patternlist:
			if complete_data[mod][gene] > 0.000001:
				f.write(mod+'\t')
		f.write('\n')
#			else:
#				f.write('\t0')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
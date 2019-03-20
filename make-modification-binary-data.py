#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
#import BED;
#import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def gene_modification_denary_set(datafile):
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene = sline[0]
			del sline[0]
			result[gene] = atoi(''.join(sline),2)
	file.close()
	return result


def gene_modification_binary_set(datafile):
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene = sline[0]
			del sline[0]
			result[gene] = ''.join(sline)
	file.close()
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-t", "--typechoice", action="store", type="int", dest="choice", metavar="<int>", help="output type, 1 for binary data, 2 for denary data")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if not opt.choice in [1,2]:
		print "output type choice error"
	else:
		if opt.choice == 1:
			data = gene_modification_binary_set(opt.infile)
		elif opt.choice == 2:
			data = gene_modification_denary_set(opt.infile)
		f = open(opt.output,'w')
		for gene in data.keys():
			f.write(gene + '\t' + str(data[gene]) + '\n')
		f.close()


if __name__ == "__main__":
	main(sys.argv)
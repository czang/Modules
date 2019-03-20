#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;




def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number


def get_formated_genes(input_file, output_file):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	#gene_list = []
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			tmp = sline[0]
			sline[0] = sline[1]
			sline[1] = tmp
			outfile.write('\t'.join(sline) + '\n')
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	get_formated_genes(opt.input_file, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

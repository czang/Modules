#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
from gene_set_manipulation import *
from associate_binary_modification_with_expression import *



"""
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""


def normalize_tag_count(input_file, number, total, output_file):
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
			if len(sline) >= number:
				sline[number-1] = str(atof(sline[number-1]) - total)
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_float_list(opt.input_file, opt.colum-1)
	Median = middle(List)
	normalize_tag_count(opt.input_file, opt.colum, Median, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

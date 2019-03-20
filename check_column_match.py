#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
import get_total_tag_counts
import gene_set_manipulation



"""
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""


def colmatch(input_file, output_file, col1=0, col2=3, col3=1, col4=4):
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
			assert sline[col1] == sline[col2]
			if sline[col3] == sline [col4]:
				sline.append(str(1))
			else:
				sline.append(str(0))
			outline = '\t'.join(sline) +'\n'
			outfile.write(outline)
	infile.close()
	outfile.close()


def total_counts(file, column):
	total = 0.0
	infile = open(file,'r')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= column:
				total += atof(sline[column-1])
	infile.close()
	return total


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	

	colmatch(opt.input_file, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

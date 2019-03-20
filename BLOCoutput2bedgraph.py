#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import numpy

## get BED module
import BED;
#import UCSC;



"""
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""


def chromname(number):
	if number in ch.keys():
		return ch[number]
	else:
		return number


def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number


def bloc2bed(input_file, chrom, output_file):
	#value_list = []
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	#gene_list = []
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			#if sline[0] == chrom:
			outfile.write(chrom + '\t' + sline[0] + '\t' + sline[1]  + '\t1\n')
				#value_list.append(atof(sline[3]))
	infile.close()
	outfile.close()
	#return value_list


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-b", "--chromosome", action="store", type="string", dest="chrom", metavar="<string>", help="chromosome")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	bloc2bed(opt.input_file, opt.chrom, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

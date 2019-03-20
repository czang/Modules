#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;



"""
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""


def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number


def get_formated_genes(input_file, chrom, start, end, output_file):
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
			if len(sline) >= max(chrom, start, end):
				outfile.write(chromname(sline[chrom - 1]) + '\t' + sline[start - 1] + '\t' + sline[end - 1]  + '\t1\n')
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
	
	#get_formated_genes(opt.input_file, opt.chrom, opt.start, opt.end, opt.output_file)
	bed = open(opt.input_file, 'r')
	bed_out = open(opt.output_file, 'w')
	plus_color = '0,0,255'
	minus_color = '255,0,0'
	bed_out.write("track name="+opt.output_file+" itemRgb=On\n")
	for line in bed:
		line = line.strip()
		sline = line.split()
        	if '+' in line:
                	bed_out.write(line + '\t' + sline[1] + '\t'+ sline[2] +'\t' + plus_color + '\n')
        	else:
                	bed_out.write(line + '\t' + sline[1] + '\t'+ sline[2] +'\t' + minus_color + '\n')
	bed_out.close()
	#bed_minus.close()



if __name__ == "__main__":
	main(sys.argv)

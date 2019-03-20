#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def normalize_tag_count(input_file, number, old, new, output_file):
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
				sline[number-1] = str(atof(sline[number-1])/ new * old)
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	parser.add_option("-t", "--totalnumber", action="store", type="float", dest="total", metavar="<float>", help="total tag number to normalized by, -1 if counted from input file automatically")
	parser.add_option("-n", "--norm", action="store", type="float", dest="norm", metavar="<float>", help="new normalization factor")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.total < 0.01:
		total = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.input_file)
	else:
		total = opt.total
	
	normalize_tag_count(opt.input_file, opt.colum, total, opt.norm, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

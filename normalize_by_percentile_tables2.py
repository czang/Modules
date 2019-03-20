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



def normalize_percentile(input_file, number, gc, denominator, output_file):
	"""
	number is the 1-based column number 
	Return a list of names
	
	"""
	#assert len(numerator) == 100
	assert len(denominator) == 100
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	i = 0
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= number:
				raw = atof(sline[number-1])
				index = int(gc[i]*100)
				norm = raw / denominator[index]
				sline[number-1] = str(round(norm,5))
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
				i += 1
	print i
	infile.close()
	outfile.close()



def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	parser.add_option("-b", "--readpercentile", action="store", type="string", dest="readtable", metavar="<file>", help="percentile table for read file")
	parser.add_option("-d", "--peakpercentile", action="store", type="string", dest="peaktable", metavar="<file>", help="percentile table for peak file")
	parser.add_option("-g", "--gcdata", action="store", type="string", dest="gcdata", metavar="<file>", help="peak gc content file")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	
	#num = gene_set_manipulation.get_float_list(opt.peaktable,1)
	den = gene_set_manipulation.get_float_list(opt.readtable,1)
	gc = gene_set_manipulation.get_float_list(opt.gcdata,8)
	normalize_percentile(opt.input_file, opt.colum, gc, den, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

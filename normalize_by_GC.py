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


def get_float_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[c]))
	file.close()
	return List


def get_gene_float_dic(infile, colum):
	'''returns a dictionary with geneIDs as keys, expression value as values. colum is the number of colum -1 where the expression data are in the file'''
	file = open(infile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			dic[sline[0]] = atof(sline[colum])
	file.close()
	return dic


def normalize_percentile_dic(input_file, number, gc, denominator, output_file):
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
				index = int(gc[sline[0]]*100)
				norm = raw / denominator[index]
				sline[number-1] = str(round(norm,5))
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
				i += 1
	#print i
	infile.close()
	outfile.close()



def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 1")
	parser.add_option("-b", "--readpercentile", action="store", type="string", dest="readtable", metavar="<file>", help="percentile table for read file")
	parser.add_option("-g", "--gcdata", action="store", type="string", dest="gcdata", metavar="<file>", help="peak gc content file")
	parser.add_option("-c", "--GCcolumn", action="store", type="int", dest="columGC", metavar="<int>", help="colum number for GC for peak GC content file, start from 1")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	
	
	den = get_float_list(opt.readtable,1)
	gc = get_gene_float_dic(opt.gcdata, opt.columGC-1)
	normalize_percentile_dic(opt.input_file, opt.colum, gc, den, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

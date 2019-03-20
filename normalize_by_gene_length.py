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
import associate_binary_modification_with_expression





def normalize_by_gene_length(input_file, number, gene_file, output_file, extension = 1000):
	"""
	
	"""
	#infile = open(input_file,'r')
	gene_length_dic = {}
	genefile = open(gene_file, 'r')
	for line in genefile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_length_dic[sline[0]] = max(atof(sline[4]) - atof(sline[3]) + 1 - extension, 1)
	genefile.close()
	
	data_dic = associate_binary_modification_with_expression.get_expression_data_dic(input_file, number)
	
	outfile = open(output_file, 'w')
	for gene in data_dic.keys():
		if gene in gene_length_dic.keys():
			value = data_dic[gene] / gene_length_dic[gene] * 1000
			outfile.write(gene + '\t' + str(value) + '\n')
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-a", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for tag counts to be normalized, start from 0")
	parser.add_option("-k", "--genefile", action="store", type="string", dest="gene_file", metavar="<file>", help="gene file in UCSC known gene format")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	
	normalize_by_gene_length(opt.input_file, opt.colum, opt.gene_file, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

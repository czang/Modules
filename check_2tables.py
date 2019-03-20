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

def store_dic_list(infile):
	outdic = {}
	for line in infile:
		line = line.strip()
		sline = line.split()
		name = sline[0]
		List = []
		for i in range(1,len(sline)):
			item = atof(sline[i])
			List.append(item)
		outdic[name] = List
	infile.close()
	return outdic


def compare2(input_file, input_file2, output_file, cut1=3.0, cut2=0.01):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	infile = open(input_file,'r')
	infile2 = open(input_file2,'r')
	outfile = open(output_file, 'w')
	h1 = infile.readline()
	h2 = infile2.readline()
	h1 = h1.strip()
	types = h1.split()
	types2 = h2.split()
	assert len(types) == len(types2)
	dic1 = store_dic_list(infile)
	dic2 = store_dic_list(infile2)
	for gene in dic1.keys():
		List1 = dic1[gene]
		List2 = dic2[gene]
		assert len(List1) == len(List2)
		for i in range(0,len(List1)):
			if List1[i] >= cut1 and List2[i] <= cut2:
				outline = gene + '\t' + types[i] + '\t' + str(List1[i]) + '\t' + str(List2[i]) + '\n'
				outfile.write(outline)
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
	parser.add_option("-j", "--inputfile2", action="store", type="string", dest="input_file2", metavar="<file>", help="original data file to be normalized")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	

	compare2(opt.input_file, opt.input_file2, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

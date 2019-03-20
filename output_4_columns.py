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


def get_formated_genes(input_file, name, value1, value2, value3, output_file):
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
			if len(sline) >= max(name, value1, value2, value3):
				outfile.write(sline[name - 1] + '\t' +sline[value1 - 1] + '\t' +sline[value2 - 1] + '\t' +sline[value3 - 1] + '\n')
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-b", "--name", action="store", type="int", dest="name", metavar="<int>", help="colum number for name")
	parser.add_option("-f", "--Value1", action="store", type="int", dest="value1", metavar="<int>", help="colum number for value1")
	parser.add_option("-g", "--Value2", action="store", type="int", dest="value2", metavar="<int>", help="colum number for value2")
	parser.add_option("-k", "--Value3", action="store", type="int", dest="value3", metavar="<int>", help="colum number for value3")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	get_formated_genes(opt.input_file, opt.name, opt.value1, opt.value2, opt.value3, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

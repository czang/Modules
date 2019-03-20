#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *

Dir = os.getcwd();


def get_float_data_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	file = open(gene_file,'r')
	gene_list = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_list.append(log(atof(sline[c]),10))
	file.close()
	return gene_list


def average(List):
	total = 0.0
	for item in List:
		total += item
	return total/len(List)


def middle(List):
	if len(List)!=0:
		if len(List)%2 == 1:
			return List[len(List)/2]
		else:
			return (List[len(List)/2-1]+List[len(List)/2])/2.0
	else: return 0.0


def standard_deviation(List):
	assert len(List)!=0
	total = 0.0
	squared_total = 0.0
	for item in List:
		total += item
		squared_total += pow(item,2)
	return sqrt(squared_total/len(List) - pow(total/len(List),2))


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name in USCS format")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	datalist = get_float_data_list(opt.infile, 12)
	print len(datalist), average(datalist), standard_deviation(datalist)

if __name__ == "__main__":
	main(sys.argv)
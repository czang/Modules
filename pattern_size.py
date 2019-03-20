#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

Dir = os.getcwd();


def get_data_list(gene_file, c):
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
			gene_list.append(sline[c])
	file.close()
	return gene_list


def histogram_data(datalist):
	datalist.sort()
	datalist.reverse()
	histogram = {}
	i = 0
	while i < (len(datalist)-1):
		j = i+1
		histogram[datalist[i]] = 1
		while j < len(datalist) and datalist[j] == datalist[i]:
			histogram[datalist[i]]+=1
			j+=1
		i = j
	return histogram


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	decimallist = get_data_list(opt.infile, 1)
	histogramlist = histogram_data(decimallist)
	f = open(opt.output,'w')
	for i in histogramlist.keys():
			f.write(str(i) + '\t' + str(histogramlist[i]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
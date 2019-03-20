#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *

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

def get_list(gene_file, c):
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


def combinatorial_pattern_gene_set(namelist, datalist, setsize):
	assert len(namelist) == len(datalist)
	#datalist.sort()
	#datalist.reverse()
	whole_set = {}
	i = 0
	while i < (len(datalist)-1):
		j = i+1
		k = 1
		while j < len(datalist) and datalist[j] == datalist[i]:
			k+=1
			j+=1
		if k >= setsize:
			subset = []
			for index in range(i, i+k):
				subset.append(namelist[index])
			whole_set[datalist[i]] = subset
		i = j
	return whole_set


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-n", "--genenumber", action="store", type="int", dest="threshold", metavar="<int>", help="minumum gene number in a combinatorial pattern")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	decimallist = get_data_list(opt.infile, 1)
	genelist = get_list(opt.infile, 0)
	whole_dic = combinatorial_pattern_gene_set(genelist, decimallist, opt.threshold)
	for pattern in whole_dic.keys():
		output_UCSCsubset_in_file (opt.known_genes, whole_dic[pattern], opt.outfile+'_'+str(pattern)+'_'+str(len(whole_dic[pattern])))


if __name__ == "__main__":
	main(sys.argv)
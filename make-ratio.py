#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;

import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is used to compare two sets of gene names and find the same and different genes between them. 
The input are two files, with the RefSeq IDs of the genes in the first column of the first file and in the second column of the second file. 
The output are three files of gene RefSeq lists, which are the same genes and the different genes for each one. 

"""


def getNameLists(IDfile, c):
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		line = line.strip()
		sline = line.split()
		refseq[sline[0]] = atof(sline[c])
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--genefile1", action="store", type="string", dest="genefile1", metavar="<file>", help="file with a list")
	parser.add_option("-b", "--genefile2", action="store", type="string", dest="genefile2", metavar="<file>", help="file with b list")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name of a/b")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.output,'w')
	data_dic_1 = getNameLists(opt.genefile1, 1)
	data_dic_2 = getNameLists(opt.genefile2, 1)
	for gene in data_dic_1.keys():
		if data_dic_2.has_key(gene):
			f.write(gene + '\t' + str(data_dic_1[gene] / data_dic_2[gene]) + '\n')


if __name__ == "__main__":
	main(sys.argv)
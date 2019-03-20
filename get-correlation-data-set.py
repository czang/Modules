#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is to combine two data files with their GeneIDs, generating a file of two columns of data so that it can be used to plot the correlation figure. GeneIDs are still in the first column.
"""


def getNameLists(IDfile, c):
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		line = line.strip()
		sline = line.split()
		refseq[sline[0]] = sline[c]
	file.close()
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--genedatafilex", action="store", type="string", dest="genefile1", metavar="<file>", help="file with gene expression list")
	parser.add_option("-y", "--genedatafiley", action="store", type="string", dest="genefile2", metavar="<file>", help="file with gene island tag list")
	parser.add_option("-a", "--datacolumx", action="store", type="int", dest="xcolum", metavar="<int>", help="x data colum number to combine, start from 0")
	parser.add_option("-b", "--datacolumy", action="store", type="int", dest="ycolum", metavar="<int>", help="y data colum number to combine, start from 0")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.output,'w')
	data_dic_1 = getNameLists(opt.genefile1, opt.xcolum)
	data_dic_2 = getNameLists(opt.genefile2, opt.ycolum)
	for gene in data_dic_1.keys():
		if data_dic_2.has_key(gene):
			f.write(gene + '\t' + data_dic_1[gene] + '\t' + data_dic_2[gene] + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
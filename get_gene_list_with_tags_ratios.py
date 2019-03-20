#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
import gene_set_manipulation;




def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--genefile1", action="store", type="string", dest="genefile1", metavar="<file>", help="file with gene expression info and RefSeq should be in the second columne")
	parser.add_option("-b", "--genefile2", action="store", type="string", dest="genefile2", metavar="<file>", help="file with island analysis on gene and RefSeq should be in the first columne")
	parser.add_option("-c", "--genefile3", action="store", type="string", dest="genefile3", metavar="<file>", help="file with island analysis on gene and RefSeq should be in the first columne")
	parser.add_option("-n", "--number", action="store", type="int", dest="number", metavar="<int>", help="column number, start from 0")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="file with island analysis on gene and RefSeq should be in the first columne")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	f1 = gene_set_manipulation.get_gene_dictionary(opt.genefile1, opt.number)
	f2 = gene_set_manipulation.get_gene_dictionary(opt.genefile2, opt.number)
	f3 = gene_set_manipulation.get_gene_dictionary(opt.genefile3, opt.number)
	
	f = open(opt.outfile,'w')
	for gene in f1.keys():
		f.write(gene + '\t' + f1[gene] + '\t' + f2[gene] + '\t' + f3[gene] + '\n')
	f.close()
	result = {}



if __name__ == "__main__":
	main(sys.argv)

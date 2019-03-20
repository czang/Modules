#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import gene_set_manipulation
import separate_by_chrom
import associate_binary_modification_with_expression

"""
Combine the tag-ratio data with the gene expression data
"""

#def 

def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--genefile1", action="store", type="string", dest="genefile1", metavar="<file>", help="file with gene expression list")
	parser.add_option("-y", "--genefile2", action="store", type="string", dest="genefile2", metavar="<file>", help="file with gene tag ratio list")
	#parser.add_option("-c", "--datacolum", action="store", type="int", dest="number", metavar="<int>", help="data colum number to combine")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	data_dic_1 = gene_set_manipulation.get_gene_dictionary(opt.genefile1, 1)
	data_dic_2 = gene_set_manipulation.get_gene_dictionary(opt.genefile2, 1)
	
	f = open(opt.outfile,'w');
	for gene in data_dic_1.keys():
		if data_dic_2.has_key(gene):
			f.write(gene + '\t' + data_dic_1[gene] + '\t' + data_dic_2[gene] + '\n')
	f.close();

if __name__ == "__main__":
	main(sys.argv)
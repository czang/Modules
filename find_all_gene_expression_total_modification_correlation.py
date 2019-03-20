#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--datafiletail", action="store", type="string", dest="infile", metavar="<file>", help="input the common tail of modification data files")
	parser.add_option("-p", "--datapath", action="store", type="string", dest="datapath", metavar="<file>", help="input the path where modification data files are stored")
	parser.add_option("-k", "--expressionfile", action="store", type="string", dest="expressionfile", metavar="<file>", help="input geneset file with expression info")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	modificationlist = get_complete_modification_level_list(opt.infile, opt.datapath)
	expressiondic = get_expression_log_data_dic(opt.expressionfile, 12)
	modification_level_dic = get_gene_total_modification_level_sum_Log(modificationlist, expressiondic.keys())
	
	f = open(opt.output, 'w')
	for gene in expressiondic.keys():
		assert modification_level_dic.has_key(gene)
		f.write(gene +'\t'+ str(expressiondic[gene]) +'\t'+ str(modification_level_dic[gene]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
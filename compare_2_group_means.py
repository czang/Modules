#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile1", metavar="<file>", help="input file with control group")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="inputfile2", metavar="<file>", help="input file with target group")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	#Average_List = []
	#Average_List_2 = []
	#data_dic = {}
	#data_dic_2 = {}
	data_dic = get_positive_log_data_dic(opt.inputfile1, 1)
	data_dic_2 = get_positive_log_data_dic(opt.inputfile2, 1)
	Average_1 = pow(10,average_expression_of_genes(data_dic, data_dic.keys()))
	Average_2 = pow(10,average_expression_of_genes(data_dic_2, data_dic_2.keys()))
	ratio = Average_2/Average_1
	
	f = open(opt.out_file, 'w')
	f.write(str(Average_1)+'\t'+str(Average_2)+'\t'+str(ratio)+'\n')
	f.close()
	print opt.out_file, Average_1, Average_2, ratio


if __name__ == "__main__":
	main(sys.argv)
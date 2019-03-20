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
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-n", "--groupnumber", action="store", type="int", dest="groupnumber", metavar="<int>", help="number of groups")
	parser.add_option("-c", "--columnnumber", action="store", type="int", dest="column", metavar="<int>", help="number of column of the expression data, start from 0")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	Average_List = []
	Deviation_List = []
	data_dic = {}
	for i in range(0, opt.groupnumber):
		data_dic = get_expression_log_data_dic(opt.inputfile+'_'+str(i+1), opt.column)
		Average_List.append(average_expression_of_genes(data_dic, data_dic.keys()))
		Deviation_List.append(expression_deviation_of_genes(data_dic, data_dic.keys()))
	
	f = open(opt.out_file, 'w')
	for i in range(0, len(Average_List)):
		f.write(str(i+1)+'\t'+str(Average_List[i])+'\t'+str(Deviation_List[i])+'\n')
	f.close()



if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

number_of_modifications = 42


def get_result_dictionary(numberdic, expressiondic):
	'''returns a list, index as the number of markers, value as the average value of expression levels'''
	List = []
	for i in range(0,number_of_modifications):
		group = []
		for gene in numberdic.keys():
			if numberdic[gene] == i and expressiondic.has_key(gene):
				group.append(expressiondic[gene])
		if len(group) != 0:
			a = 0.0
			for item in group:
				a += item
			List.append(a / len(group))
			print i, List[i], len(group)
		else:
			List.append(0)
	return List


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--completedata", action="store", type="string", dest="infile", metavar="<file>", help="input whole data set file name")
	parser.add_option("-k", "--genefile", action="store", type="string", dest="expressionfile", metavar="<file>", help="input gene set file with expression info")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	datadic = read_complete_binary_file(opt.infile)
	expressiondic = get_expression_log_data_dic(opt.expressionfile, 12)
	markernumberdic = get_modification_number_for_genes(datadic, expressiondic.keys())
	result = get_result_dictionary(markernumberdic, expressiondic)
	
	f = open(opt.output, 'w')
	for i in range(0, len(result)):
		if result[i] != 0:
			f.write(str(i)+'\t'+str(result[i])+'\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
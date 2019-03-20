#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--file1", action="store", type="string", dest="file1", metavar="<file>", help="file 1")
	parser.add_option("-b", "--file2", action="store", type="string", dest="file2", metavar="<file>", help="file 2")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	dic1 = get_data_dic(opt.file1, 0, 1)
	dic2 = get_data_dic(opt.file2, 0, 1)
	comparison = gene_comparison(dic1.keys(), dic2.keys())
	f = open(opt.out_file+'_shared', 'w')
	for item in comparison["shared"]:
		f.write(item + '\n')
	f.close()
	f = open(opt.out_file+'_only1', 'w')
	for item in comparison['only in 1']:
		f.write(item + '\t' + str(dic1[item]) + '\n')
	f.close()
	f = open(opt.out_file+'_only2', 'w')
	for item in comparison['only in 2']:
		f.write(item + '\t' + str(dic2[item]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
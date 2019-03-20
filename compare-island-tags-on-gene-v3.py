#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
import separate_by_chrom
from associate_binary_modification_with_expression import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def dic_log_change(dic_1, dic_2, logbase = 2):
	result = {}
	if len(dic_1.keys()) != len(dic_2.keys()):
		print 'list lengths not equal!'
	else:
		for name in dic_1.keys():
			if dic_2.has_key(name):
				if dic_1[name] < 1e-7:
					result[name] = 0.0
					if dic_2[name] > 1e-7:
						result[name] = result[name] + log(dic_2[name],logbase)
				else:
					result[name] = -log(dic_1[name],logbase)
					if dic_2[name] > 1e-7:
						result[name] = result[name] + log(dic_2[name],logbase)
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile1", action="store", type="string", dest="file1", metavar="<file>", help="initial state file to be denominator")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="file2", metavar="<file>", help="final state file to be numerator")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	dic_1 = get_expression_data_dic(opt.file1, 1)
	dic_2 = get_expression_data_dic(opt.file2, 1)
	Result = dic_log_change(dic_1, dic_2)
	
	f = open (opt.outfile,'w')
	for g in Result.keys():
		f.write(g + '\t' + str(Result[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
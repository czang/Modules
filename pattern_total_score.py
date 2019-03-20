#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import associate_binary_modification_with_expression
from gene_set_manipulation import *

cat = "/bin/cat"


fib_index_list = [0, 1, 2, 3, 4, 5, 12, 13, 18, 31, 33, 36]
epi_index_list = [15, 16, 21, 23, 24, 25, 26, 32, 34, 35]
blood_index_list = [7, 8, 9, 20, 27, 28, 30, 38, 39]

def index_score(pattern_list, index_list):
	score = 0
	for i in index_list:
		assert i < len(pattern_list)
		score += pattern_list[i]
	return score


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="pattern file name")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	pattern_dic = associate_binary_modification_with_expression.read_complete_binary_file(opt.infile)
	f = open(opt.infile, 'r')
	o = open(opt.outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[0]
			List = []
			for i in range(1, len(sline)):
				List.append(atoi(sline[i]))
			#blood_score = index_score(List, blood_index_list)
			#epi_score = index_score(List, epi_index_list)
			#fib_score = index_score(List, fib_index_list)
			total_score = index_score(List, range(0,len(sline)-1))
			o.write(name + '\t' + str(total_score) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
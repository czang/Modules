#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

## get BED module
import BED
from GenomeData import *
from associate_binary_modification_with_expression import *
#import separate_by_chrom

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");


def correlation1(list1, list2):
	assert len(list1) == len(list2)
	N = len(list1)
	avg1 = average(list1)
	avg2 = average(list2)
	total = 0.0
	for i in range(0,N):
		total += (list1[i] - avg1) * (list2[i] - avg2)
	total = total / N
	return total


def get_data_list(file, colum):
	'''colum starts from 0'''
	file = open(file, 'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[colum]))
	file.close()
	return List


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--data file", action="store", type="string", dest="datafile", metavar="<file>", help="data file")
	#parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	
	List1 = get_data_list(opt.datafile, 1)
	List2 = get_data_list(opt.datafile, 2)
	result = correlation1(List1, List2) / sqrt(correlation1(List1, List1)) / sqrt(correlation1(List2, List2))
	print opt.datafile, result


if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
import stats

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


def group_list_geo(List, size):
	result = []
	length = len(List)
	i = 0
	while i < length:
		templist = []
		j = 0
		while j < size and i < length:
			templist.append(log(List[i]))
			j += 1
			i += 1
		result.append(stats.mean(templist))
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--data file", action="store", type="string", dest="datafile", metavar="<file>", help="data file")
	parser.add_option("-n", "--groupsize", action="store", type="int", dest="groupsize", metavar="<int>", help="number of elements in each group")

	#parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	
	rawList1 = get_data_list(opt.datafile, 1)
	rawList2 = get_data_list(opt.datafile, 2)
	List1 = group_list_geo(rawList1, opt.groupsize)
	List2 = group_list_geo(rawList2, opt.groupsize)
	#result = correlation1(List1, List2) / sqrt(correlation1(List1, List1)) / sqrt(correlation1(List2, List2))
	result = stats.lpearsonr(List1, List2)[0]
	print opt.datafile, result


if __name__ == "__main__":
	main(sys.argv)

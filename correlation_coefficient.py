#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import numpy


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def get_correlation_list(datafile, column):
	file = open(datafile,'r')
	result = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result.append(atof(sline[column]))
	file.close()
	return result


def average(List):
	total = 0.0
	for item in List:
		total += item
	return total/len(List)


def middle(List):
	if len(List)!=0:
		if len(List)%2 == 1:
			return List[len(List)/2]
		else:
			return (List[len(List)/2-1]+List[len(List)/2])/2.0
	else: return 0.0


def standard_deviation(List):
	total = 0.0
	squared_total = 0.0
	for item in List:
		total += item
		squared_total += pow(item,2)
	return sqrt(squared_total/len(List) - pow(total/len(List),2))


def merge_lists(List1, List2):
	List = List1
	for item in List2:
		if not item in List1:
			List.append(item)
	return List


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input_file", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-a", "--column_number1", action="store", type="int", dest="cla", metavar="<int>", help="input column number of the first data list, start from 0")
	parser.add_option("-b", "--olumn_number2", action="store", type="int", dest="clb", metavar="<int>", help="input column number of the second data list, start from 0")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	list1 = get_correlation_list(opt.infile, opt.cla)
	list2 = get_correlation_list(opt.infile, opt.clb)
	print numpy.corrcoef(list1,list2)[0][1]
	'''f = open(opt.output,'w')
	for mod in modifications:
		f.write(mod + '\t' + str(get_mean_expression_change_with_modification(datadic, expressiondic, mod)) + '\n')
	f.close()'''


if __name__ == "__main__":
	main(sys.argv)
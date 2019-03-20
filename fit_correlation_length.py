#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
import stats
import numpy
import scipy.stats

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


def correlation_length_fit(xlist, ylist):
	assert len(xlist) == len(ylist)
	loglist = []
	s = ylist[0] * ylist[0]
	for i in range(0, len(ylist)):
		loglist.append(log(max(ylist[i] - s, 0.000000000001)))
	(a,b,r,stderr,x) = scipy.stats.linregress(xlist,loglist)
	return -1.0/a


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
	parser.add_option("-x", "--xcolumn", action="store", type="int", dest="x", help="x colomn, starting with 0", metavar="<int>")
	parser.add_option("-y", "--ycolumn", action="store", type="int", dest="y", help="y colomn, starting with 0", metavar="<int>")
	parser.add_option("-n", "--datapoints", action="store", type="int", dest="n", help="number of data points, >= 2", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	
	List1 = get_data_list(opt.datafile, opt.x)
	List2 = get_data_list(opt.datafile, opt.y)
	assert len(List1) == len(List2)
	if opt.n < len(List1):
		xlist = List1[0:(opt.n - 1)]
		ylist = List2[0:(opt.n - 1)]
	else:
		xlist = List1
		ylist = List2
	#result = correlation1(List1, List2) / sqrt(correlation1(List1, List1)) / sqrt(correlation1(List2, List2))
	#result2 = stats.lpearsonr(List1, List2)[0]
	print opt.datafile, correlation_length_fit(xlist, ylist)


if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from stats import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
"""


def get_index(List, x):
	List.sort()
	y = len(List)
	i = 0
	while i < len(List):
		if x <= List[i]:
			y = i
			break
		else:
			i += 1
	return y


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="infile", metavar="<file>", help="file with data set")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)

	f = open(opt.infile, 'r')
	o = open(opt.output, 'w')
	for line in f:
		line = line.strip()
		sline = line.split()
		name = sline[0]
		List = []
		for i in range(1, len(sline)):
			List.append(atof(sline[i]))
		flag = median(List) + 2.0
		SortList = List
		SortList.sort()
		ind = get_index(SortList, flag)
		if ind < len(List) and ind >= len(List) - 4 and (flag - SortList[0]) <= 3.0: 
			List = []
			for i in range(1, len(sline)):
				List.append(atof(sline[i]))
			index_list = []
			for i in range(ind, len(List)):
				score = SortList[i]
				index = List.index(score)
				index_list.append(str(index+1))
			o.write(name + '\t' + '\t'.join(index_list) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
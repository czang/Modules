#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import stats


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to calculate the pearson correlation between gene expression and a group of DNase peaks. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--sample1", action="store", type="string", dest="expressionfile", metavar="<file>", help="expression file name")
	parser.add_option("-b", "--sample2", action="store", type="string", dest="peakfile", metavar="<file>", help="Peak file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	"Has to be in one gene!!!"
	

	t = open(opt.expressionfile, 'r')
	expression_list = []
	l = t.readline()
	l = l.strip()
	sl = l.split()
	for i in range(1,len(sl)):
		expression_list.append(atof(sl[i]))
	t.close()
	
	p = open(opt.peakfile, 'r')
	u = open(opt.out_file, 'w')
	for line in p:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			peakList = []
			for i in range(4,len(sline)):
				peakList.append(atof(sline[i]))
			assert len(peakList) == len(expression_list)
			r = stats.lpearsonr(peakList, expression_list)[0]
			u.write(sline[0] + '\t' + str(r) + '\n')
	p.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
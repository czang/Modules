#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from stats import *
import bisect
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--motif", action="store", type="string", dest="motif", metavar="<file>", help="motif ID")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	List = []
	inputfile = open(opt.motif+'.txt','r')
	for line in inputfile:
		if not (re.match("#", line) or re.match("track", line) or re.match("variableStep", line)):
			line = line.strip()
			sline = line.split()
			assert len(sline) >= 5
			List.append(atof(sline[4]))
	inputfile.close()
	Length = len(List)
	List.sort()
	nonZero = Length - bisect.bisect_left(List, 0.5)
	if nonZero > 0: 
		Sum = sum(List)
		n0Mean = sum(List) / nonZero
		n0Median = List[int(Length - nonZero / 2)]
		print opt.motif, nonZero, Sum, n0Mean, n0Median
	else:
		print opt.motif, 0, 0, 0, 0
	
if __name__ == "__main__":
	main(sys.argv)

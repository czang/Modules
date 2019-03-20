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


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="infile", metavar="<file>", help="file with data set")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)

	f = open(opt.infile, 'r')
	dict = {}
	for line in f:
		line = line.strip()
		sline = line.split()
		name = sline[0]
		for i in range(1, len(sline)):
			cellid = sline[i]
			if cellid in dict.keys():
				dict[cellid].append(name)
			else:
				List = []
				List.append(name)
				dict[cellid] = List
	f.close()
	for cellid in dict.keys():
		o = open(opt.output+cellid, 'w')
		List = dict[cellid]
		for item in List:
			o.write(item + '\n')
		o.close()


if __name__ == "__main__":
	main(sys.argv)
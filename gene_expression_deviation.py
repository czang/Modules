#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from associate_binary_modification_with_expression import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="infile", metavar="<file>", help="sample name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.infile, 'r')
	u = open(opt.out_file, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			datalist = []
			for i in range(1,len(sline)):
				datalist.append(atof(sline[i]))
			dev = standard_deviation(datalist)
			u.write(sline[0] + '\t' + str(dev) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
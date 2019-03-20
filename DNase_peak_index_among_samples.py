#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	parser.add_option("-c", "--sample", action="store", type="int", dest="index", metavar="<int>", help="sample label")
	parser.add_option("-i", "--file1", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.infile, 'r')
	o = open(opt.outfile, 'w')
	for line in f:
		if not re.match("#",line):
			line = line.strip()
			sline = line.split()
			assert opt.index <= len(sline)
			expression = atof(sline[opt.index + 3])
			List = []
			for i in range(4,len(sline)):
				List.append(atof(sline[i]))
			List.sort()
			index = List.index(expression) + 1
			o.write(sline[0] + '\t' + str(index) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--file1", action="store", type="string", dest="file1", metavar="<file>", help="file 1")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.file1, 'r')
	o = open(opt.outfile, 'w')
	for line in f:
		if not re.search("A", line):
			line = line.strip()
			sline = line.split()
			for i in range(0, len(sline)):
				sline[i] = str(atof(sline[i]) * 100)
			o.write('\t'.join(sline) + '\n')


if __name__ == "__main__":
	main(sys.argv)
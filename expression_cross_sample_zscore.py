#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from stats import *


def main(argv):
	parser = OptionParser()
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
			List = []
			for i in range(1,len(sline)):
				List.append(atof(sline[i]))
			mu = mean(List)
			sigma = stdev(List)
			Result = []
			for i in range(0, len(List)):
				Result.append(str((List[i] - mu)/sigma))
			o.write(sline[0] + '\t' + '\t'.join(Result) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
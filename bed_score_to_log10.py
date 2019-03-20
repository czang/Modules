#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="inputfile", metavar="<file>", help="input BED file")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.inputfile, 'r')
	o = open(opt.outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			assert len(sline) > 5
			score = log(atof(sline[4]),10)
			sline[4] = str(round(score,4))
			o.write('\t'.join(sline) + '\n')
	f.close()
	o.close()

if __name__ == "__main__":
	main(sys.argv)       
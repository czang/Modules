#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


def remove_lines(infile, outfile, c = 7):
	f = open(infile,'r')
	o = open(outfile,'w')
	i = 0
	for line in f:
		line = line.strip()
		sline = line.split()
		if len(sline) == 7: 
			o.write('\t'.join(sline) + '\n')
		else:
			i += 1
	f.close()
	o.close()
	return i


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	print remove_lines(opt.infile, opt.outfile)


if __name__ == "__main__":
	main(sys.argv)
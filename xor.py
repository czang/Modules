#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator


def xor2(infile, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			x = atoi(sline[1])
			y = atoi(sline[2])
			z = x ^ y
			sline[2] = str((z-x+y)/2)
			o.write('\t'.join(sline)+'\n')
	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--in_file", action="store", type="string",
                      dest="in_file", help="input file", metavar="<file>")
	parser.add_option("-o", "--output_file", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	xor2(opt.in_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

Dir = os.getcwd();

'''Input files must be sorted'''

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.infile, 'r')
	o = open(opt.output, 'w')
	l = f.readline()
	l = l.strip()
	sl = l.split()
	x = sl[1]
	y = sl[2]
	n = 1
	for line in f:
		line = line.strip()
		sline = line.split()
		if sline[1] == x and sline[2] == y:
			n += 1
		else:
			o.write(x + '\t' + y + '\t' + str(n) + '\n')
			x = sline[1]
			y = sline[2]
			n = 1
	o.write(x + '\t' + y + '\t' + str(n) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
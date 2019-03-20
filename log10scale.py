#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	o = open(opt.output,'w')
	f = open(opt.xfile, 'r')
	for line in f:
		line = line.strip()
		sline = line.split()
		linear = int(sline[1])
		if linear > 0:
			sline[1] = str(log(linear,10))
			o.write('\t'.join(sline) + '\n')
		elif linear < 0:
			sline[1] = str(-log(-linear,10))
			o.write('\t'.join(sline) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
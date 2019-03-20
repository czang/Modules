#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
import gene_set_manipulation;


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="List of island coordinate file")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output BED format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	window = 1000
	start = 1000
	file = open(opt.input_file,'r')
	o = open(opt.output_file, 'w')
	l = 0
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			l +=1
			x = start
			for item in sline:
				if item == '1':
					o.write('chr' + str(l) + '\t' + str(x) + '\t' + str(x+window-1) + '\t1\n')
				x += window
	file.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)

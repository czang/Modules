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
	parser.add_option("-b", "--inputfile2", action="store", type="string", dest="input_file2", metavar="<file>", help="List of island size file")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output BED format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	colum = 1
	mer = 25
	List=[]
	ref = open(opt.input_file2, 'r')
	for line in ref:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline) > colum):
				List.append(int(float(sline[colum])))
	ref.close()	
	
	file = open(opt.input_file,'r')
	o = open(opt.output_file, 'w')
	i = 0
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline) > colum):
				item = int(float(sline[colum]))
				length = List[i]
				o.write('chr1\t' + str(item) + '\t' + str(item + length) + '\t1\n')
				i += 1
	file.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)

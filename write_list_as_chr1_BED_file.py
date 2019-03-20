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
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	colum = 1
	mer = 25
	List=[]
	file = open(opt.input_file,'r')
	o = open(opt.output_file, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline) > colum):
				if not 'e' in sline[colum]:
					item = int(sline[colum])
					o.write('chr1\t' + str(item) + '\t' + str(int(item) + mer - 1) + '\t0\t0\t+\n')
	file.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)

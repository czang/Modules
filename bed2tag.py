#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;

plus = re.compile('\+');
minus = re.compile('\-');


"""
This module is to convert bed file to tag file. 

"""


def BED_to_tag(input_file, output_file):
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	#gene_list = []
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if plus.match(sline[5]):
				outfile.write(sline[0]+'\t'+sline[1]+'\t'+sline[5]+'\n')
			elif minus.match(sline[5]):
				outfile.write(sline[0]+'\t'+sline[2]+'\t'+sline[5]+'\n')
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	BED_to_tag(opt.input_file, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

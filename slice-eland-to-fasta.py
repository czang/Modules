#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

#from gene_set_manipulation import *
#import BED

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def slice_eland_to_fasta(infile, outfile, length):
	file = open(infile,'r')
	output = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			output.write(sline[0] + '\n')
			output.write(sline[1][:length] + '\n')
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-l", "--sequencelength", action="store", type="int", dest="length", help="length of output sequences", metavar="<int>")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	slice_eland_to_fasta(opt.inputfile, opt.out_file, opt.length)


if __name__ == "__main__":
	main(sys.argv)
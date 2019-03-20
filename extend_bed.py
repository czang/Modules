#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def summit_region_to_bed(infile, outfile, extension = 50):
	file = open(infile,'r')
	output = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline)>0):
				chrom = sline[0]
				start = atoi(sline[1]) - extension
				end = atoi(sline[2]) + extension
				sline[1] = str(start)
				sline[2] = str(end)
				output.write('\t'.join(sline) + '\n')
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-w", "--extension", action="store", type="int", dest="extension", help="extension", metavar="<int>")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	summit_region_to_bed(opt.inputfile, opt.out_file, opt.extension)


if __name__ == "__main__":
	main(sys.argv)
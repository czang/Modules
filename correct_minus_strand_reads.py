#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

#import BED

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def correct_minus(infile, outfile):
	file = open(infile,'r')
	output = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			assert len(sline) == 6
			if plus.match(sline[5]):
				output.write('\t'.join(sline) + '\n')
			elif minus.match(sline[5]):
				start = atoi(sline[1])
				end = atoi(sline[2])
				sline[1] = str(end)
				sline[2] = str(end * 2 - start)
				output.write('\t'.join(sline) + '\n')
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	correct_minus(opt.inputfile, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
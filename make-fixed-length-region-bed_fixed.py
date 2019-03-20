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


def summit_region_to_bed(infile, outfile, column = -1, width = 100):
	file = open(infile,'r')
	output = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline)>0):
				chrom = sline[0]
				start = atoi(sline[1])
				end = atoi(sline[2])
				length = end - start
				halfwidth = width / 2
				if column > 0:
					mid = start + atoi(sline[column - 1]) - 1
				else:
					mid = (start + end) / 2
				output_start = mid - halfwidth
				output_end = mid + halfwidth
				output.write(sline[0] + '\t' + str(output_start) + '\t' + str(output_end) + '\n')
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-w", "--maxwidth", action="store", type="int", dest="width", help="maximum width of output sites", metavar="<int>")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", help="column number of summits, start from 1, negative if using the middle position", metavar="<int>")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	summit_region_to_bed(opt.inputfile, opt.out_file, opt.column, opt.width)


if __name__ == "__main__":
	main(sys.argv)
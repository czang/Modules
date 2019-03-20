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


def summit_region_to_bed(infile, outfile):
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
				residue = length%120
				if length <= 230: 
					newlength = 240
				elif residue > 110:
					newlength = length - residue + 120
				else:
					newlength = length - residue
				halfwidth = newlength / 2
				mid = (start + end) / 2
				sline[1] = str(mid - halfwidth)
				sline[2] = str(mid + halfwidth - 1)
				output.write('\t'.join(sline) + '\n')
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	summit_region_to_bed(opt.inputfile, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
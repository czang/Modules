#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator


def add_tag_density(infile, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			corr = atof(sline[2])
			abso = abs(corr)
			sline.append(str(abso))
			o.write('\t'.join(sline)+'\n')
	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--island_bed_file", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	add_tag_density(opt.bed_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)

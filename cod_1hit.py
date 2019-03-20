#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;



"""


"""


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	infile = open(opt.input_file, 'r')
	outfile = open(opt.output_file, 'w')
	current = infile.readline()
	while re.match("#", current):
		current = infile.readline()
	current = current.strip()
	current_id = current.split()[0]
	current_score = atof(current.split()[5])
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			id = sline[0]
			score = atof(sline[5])
			if id != current_id:
				outfile.write(current+'\n')
				current = line
				current_id = id
				current_score = score
			elif score > current_score:
				current = line
				current_id = id
				current_score = score
	infile.close()
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)

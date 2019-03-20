#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	infile = open(opt.input_file, 'r');
	outfile = open(opt.output_file, 'w');
	for line in infile:
		if not re.match ("#rank",line):
			line = line.strip();
			sline = line.split();
			outline = sline[1] + "\t" + sline[2] + "\t" + sline[3] +  "\t" + sline[10] +  "\n";
			outfile.write(outline);

	infile.close()	
	outfile.close();

if __name__ == "__main__":
	main(sys.argv)



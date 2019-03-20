#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--wiggle", action="store", type="string", dest="wiggle_file", metavar="<file>", help="input wiggle file")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output wiggle file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.wiggle_file,'r')
	outfile = open(opt.out_file, 'w')
	for line in inputfile:
		if not (re.match("#", line) or re.match("track", line) or re.match("variableStep", line)):
			line = line.strip()
			sline = line.split()
			assert len(sline) >= 2
			outfile.write(sline[1] + '\n')
	outfile.close()
	inputfile.close()
	
if __name__ == "__main__":
	main(sys.argv)

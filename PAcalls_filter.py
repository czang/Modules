#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--data", action="store", type="string", dest="infile", metavar="<file>", help="input data file")
	parser.add_option("-m", "--min", action="store", type="int", dest="minP", metavar="<int>", help="minimum P")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
	inputfile = open(opt.infile,'r')
	outfile = open(opt.out_file, 'w')
	for line in inputfile:
		if not (re.match("#", line) or re.match("CEL", line)):
			line = line.strip()
			sline = line.split()
			assert len(sline) >= 2
			Counts = [0] * 3
			for i in range(1,len(sline)):
				if sline[i] == 'A':
					Counts[0] += 1
				elif sline[i] == 'P':
					Counts[2] += 1
				else:
					Counts[1] += 1
			if Counts[2] >= opt.minP:
				outfile.write('\t'.join(sline) + '\n')
	outfile.close()
	inputfile.close()
	
if __name__ == "__main__":
	main(sys.argv)

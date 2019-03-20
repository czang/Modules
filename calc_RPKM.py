#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import separate_by_chrom
import stats


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
Input file should be 
id	chrom	start	end	tags1(perM)	tags2(perM)	... 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--sample", action="store", type="string", dest="peakfile", metavar="<file>", help="DNase count file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	p = open(opt.peakfile, 'r')
	o = open(opt.out_file, 'w')
	for line in p:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			length = atof(sline[3]) - atof(sline[2])
			for i in range(4,len(sline)):
				sline[i] = str(atof(sline[i])/length*1000)
			o.write('\t'.join(sline)+'\n')
	p.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
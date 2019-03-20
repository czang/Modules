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
		line = line.strip()
		sline = line.split()
		if not re.search("NNN", sline[3]):
			o.write('\t'.join(sline)+'\n')
	p.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
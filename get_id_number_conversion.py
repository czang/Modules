#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import gene_set_manipulation


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to separate BED file into 2, according to the size of each region. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--prefix", action="store", type="string", dest="prefix", metavar="<file>", help="correlation data file")
	parser.add_option("-l", "--peaklist", action="store", type="string", dest="peakfile", metavar="<file>", help="correlation data file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.prefix+opt.peakfile, 'r')
	line = f.readline()
	line = line.strip()
	sline = line.split()
	data = sline[0]
	sdata = data.split('_')
	print sdata[0]+'_'+sdata[1], opt.peakfile

if __name__ == "__main__":
	main(sys.argv)
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
	parser.add_option("-d", "--DNase", action="store", type="string", dest="dnasefile", metavar="<file>", help="DNae data file")
	parser.add_option("-e", "--expression", action="store", type="string", dest="exprfile", metavar="<file>", help="expression data file")
	parser.add_option("-l", "--peaklist", action="store", type="string", dest="peakfile", metavar="<file>", help="correlation data file")
	parser.add_option("-m", "--matrix", action="store", type="string", dest="matrix_file", metavar="<file>", help="output DNase matrix")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output expression list")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	
	f = open(opt.peakfile, 'r')
	idlist = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			data = sline[0]
			sdata = data.split('_')
			idlist.append(int(sdata[2]))
			gene = sdata[0]+'_'+sdata[1]
	f.close()
	idlist.sort()
	peaklist = []
	for item in idlist:
		peaklist.append(str(item))
	
	gene_set_manipulation.output_UCSCsubset_in_file (opt.dnasefile, peaklist, opt.matrix_file)
	gene_set_manipulation.output_UCSCsubset_in_file (opt.exprfile, [gene], opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
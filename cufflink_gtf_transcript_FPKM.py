#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is to combine two data files with their GeneIDs, generating a file of two columns of data so that it can be used to plot the correlation figure. GeneIDs are still in the first column.
"""


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	o = open(opt.output,'w')
	f = open(opt.xfile, 'r')
	for line in f:
		line = line.strip()
		sline = line.split()
		if sline[2] == 'transcript': 
			assert sline[10] == 'transcript_id'
			assert sline[12] == 'FPKM'
			o.write(sline[11] + '\t' + sline[13] + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
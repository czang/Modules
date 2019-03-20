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


def getIntList(IDfile, c):
	file = open(IDfile,'r')
	List = []
	for line in file:
		line = line.strip()
		sline = line.split()
		List.append(atoi(sline[c]))
	file.close()
	return List


def clean_list(List1, List2, outfile):
	o = open(outfile, 'w')
	assert len(List1) == len(List2)
	i = 0
	current = List1[0]
	current2 = List2[0]
	j = 1
	while j < len(List1):
		if List1[j] == current:
			if List2[j] < current2:
				current2 = List2[j]
		else:
			o.write(str(current) + '\t' + str(current2) + '\n')
			current = List1[j]
			current2 = List2[j]
		j += 1
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="file with data set")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	List1 = getIntList(opt.infile, 0)
	List2 = getIntList(opt.infile, 1)
	clean_list(List1, List2, opt.outfile)


if __name__ == "__main__":
	main(sys.argv)
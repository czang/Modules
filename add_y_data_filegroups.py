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


def getNameLists(IDfile, c):
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		line = line.strip()
		sline = line.split()
		refseq[sline[0]] = sline[c]
	file.close()
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-y", "--y", action="store", type="string", dest="yfile", metavar="<file>", help="file with y data to be added")
	parser.add_option("-s", "--s", action="store", type="int", dest="start", metavar="<int>", help="start suffix")
	parser.add_option("-e", "--e", action="store", type="int", dest="end", metavar="<int>", help="end suffix")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	for i in range(opt.start,opt.end):
		j = str(i)
		if len(j) == 1:
			j = '00'+j
		elif len(j) == 2:
			j = '0'+j
		print j
		xfile = opt.xfile+j
		yfile = opt.yfile+j
		output = opt.output+j
		data_dic = getNameLists(yfile, 1)
		o = open(output,'w')
		f = open(xfile, 'r')
		for line in f:
			line = line.strip()
			sline = line.split()
			sline.append(data_dic[sline[0]])
			o.write('\t'.join(sline) + '\n')
		f.close()
		o.close()


if __name__ == "__main__":
	main(sys.argv)
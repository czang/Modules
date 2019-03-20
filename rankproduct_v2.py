#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is to combine two data files with their GeneIDs, generating a file of two columns of data so that it can be used to plot the correlation figure. GeneIDs are still in the first column.
"""


def getRankDic_2(IDfile):
	L = []
	file = open(IDfile,'r')
	for line in file:
		line = line.strip()
		sline = line.split()
		L.append(-1 * atof(sline[1]))
	L.sort()
	file.close()
	refseq = {}
	file = open(IDfile,'r')
	for line in file:
		line = line.strip()
		sline = line.split()
		score = -1 * atof(sline[1])
		i_f = bisect.bisect_left(L, score) + 1
		i_r = bisect.bisect_right(L, score) 
		i = float(i_f + i_r) / 2
		refseq[sline[0]] = i
	file.close()
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-y", "--y", action="store", type="string", dest="yfile", metavar="<file>", help="file with y data to be added")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	dic1 = getRankDic_2(opt.xfile)
	dic2 = getRankDic_2(opt.yfile)
	o = open(opt.output,'w')
	for gene in dic1.keys():
		if gene in dic2.keys():
			score = sqrt(float(dic1[gene] * dic2[gene]))
			o.write(gene + '\t' + str(score) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)
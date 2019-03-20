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


def getRankDic_2(IDfile,column,dr):
	L = []
	file = open(IDfile,'r')
	for line in file:
		line = line.strip()
		sline = line.split()
		L.append(dr * atof(sline[column]))
	L.sort()
	file.close()
	refseq = {}
	file = open(IDfile,'r')
	for line in file:
		line = line.strip()
		sline = line.split()
		score = dr * atof(sline[column])
		i_f = bisect.bisect_left(L, score) + 1
		i_r = bisect.bisect_right(L, score) 
		i = float(i_f + i_r) / 2
		refseq[sline[0]] = i
	file.close()
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-c", "--tagcountcolumn", action="store", type="int", dest="colum", metavar="<int>", help="colum number for data, start from 0", default=1)
	parser.add_option("-s", "--direction", action="store", type="int", dest="direction", metavar="<int>", help="1 for increasing, -1 for decreasing", default=1)
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	dic = getRankDic_2(opt.xfile,opt.colum,opt.direction)
	o = open(opt.output,'w')
	file = open(opt.xfile,'r')
	for line in file:
                line = line.strip()
                sline = line.split()
		r = dic[sline[0]]
		sline.append(str(r))
		o.write('\t'.join(sline)+'\n')
	o.close()
	file.close()


if __name__ == "__main__":
	main(sys.argv)

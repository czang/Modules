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


def getRankDic(IDfile):
	i = 1
	n = 0
	file = open(IDfile,'r')
	refseq = {}
	f = file.readline()
	f = f.strip()
	sf = f.split()
	refseq[sf[0]] = i
	current = atof(sf[1])
	n += 1
	for line in file:
		line = line.strip()
		sline = line.split()
		score = atof(sline[1])
		if score < current:
			i = n + 1
			refseq[sline[0]] = i
			current = score
			n += 1
		elif score == current:
			refseq[sline[0]] = i
			n += 1 
	file.close()
	return refseq


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--x", action="store", type="string", dest="xfile", metavar="<file>", help="file with data set")
	parser.add_option("-y", "--y", action="store", type="string", dest="yfile", metavar="<file>", help="file with y data to be added")
	parser.add_option("-z", "--z", action="store", type="string", dest="zfile", metavar="<file>", help="file with z data to be added")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	dic1 = getRankDic(opt.xfile)
	dic2 = getRankDic(opt.yfile)
	dic3 = getRankDic(opt.zfile)
	o = open(opt.output,'w')
	for gene in dic1.keys():
		if gene in dic2.keys() and gene in dic3.keys():
			score = pow((float(dic1[gene] * dic2[gene] * dic3[gene])),float(1.0/3))
			o.write(gene + '\t' + str(score) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)
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


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to separate BED file into 2, according to the size of each region. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--sample1", action="store", type="string", dest="TSSfile", metavar="<file>", help="TSS file name")
	parser.add_option("-b", "--sample2", action="store", type="string", dest="peakfile", metavar="<file>", help="Peak file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file dir")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	"Has to be in one chromosome!!!"
	
	
	p = open(opt.peakfile, 'r')
	summit_dic = {}
	for line in p:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			summit_dic[sline[0]] = (int(sline[2]) + int(sline[3]))/2
	p.close()
	
	t = open(opt.TSSfile, 'r')
	TSS_dic = {}
	for line in t:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			TSS_dic[sline[0]] = int(sline[2])
	t.close()

	for gene in TSS_dic.keys():
		u = open(opt.out_file+"/"+gene+"_dis.txt", 'w')
		for peak in summit_dic.keys():
			distance = abs(summit_dic[peak]-TSS_dic[gene])
			u.write(peak+'\t'+str(distance)+'\n')
		u.close()
	

if __name__ == "__main__":
	main(sys.argv)
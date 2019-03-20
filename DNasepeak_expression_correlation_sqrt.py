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


'''
This module is to separate BED file into 2, according to the size of each region. 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--sample1", action="store", type="string", dest="TSSfile", metavar="<file>", help="TSS file name")
	parser.add_option("-b", "--sample2", action="store", type="string", dest="peakfile", metavar="<file>", help="Peak file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	"Has to be in one chromosome!!!"
	
	
	p = open(opt.peakfile, 'r')
	DNase_dic = {}
	for line in p:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			List = []
			for i in range(5,len(sline)):
				List.append(sqrt(atof(sline[i])/atof(sline[4])*1000))
			DNase_dic[sline[0]] = List
	p.close()
	
	t = open(opt.TSSfile, 'r')
	expr_dic = {}
	for line in t:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			List = []
			for i in range(1,len(sline)):
				List.append(sqrt(atof(sline[i])))
			expr_dic[sline[0]] = List
	t.close()

	u = open(opt.out_file, 'w')
	for gene in expr_dic.keys():
		for peak in DNase_dic.keys():
			r = stats.kendalltau(DNase_dic[peak], expr_dic[gene])[0]
			u.write(gene+'_'+peak+'\t'+str(r)+'\n')
	u.close()


if __name__ == "__main__":
	main(sys.argv)
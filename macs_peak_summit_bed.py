#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to ...

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="bedfile", metavar="<file>", help="sample name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	p = open(opt.bedfile+"_summits.bed", 'r')
	summit_dic = {}
	for line in p:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			summit_dic[sline[3]] = int(sline[2])
	p.close()
	
	f = open(opt.bedfile+"_peaks.bed", 'r')
	u = open(opt.out_file, 'w')
	for line in f:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			center = summit_dic[sline[3]]
			start = center - 1
			u.write('\t'.join(sline) + '\t+\t' + str(start) + '\t' + str(center) + '\n')
	f.close()
	u.close()


if __name__ == "__main__":
	main(sys.argv)
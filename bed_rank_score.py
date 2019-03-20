#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


def organize_bed(infile, name, total, outfile):
	f = open(infile,'r')
	o = open(outfile,'w')
	i = 1
	for line in f:
		line = line.strip()
		sline = line.split()
		sline[3] = name
		sline[4] = str(float(i)/total)
		o.write('\t'.join(sline) + '\n')
		i+=1
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-n", "--name", action="store", type="string", dest="name", metavar="<string>", help="sample name")
	parser.add_option("-c", "--count", action="store", type="int", dest="total", help="total line count", metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	organize_bed(opt.infile, opt.name, opt.total, opt.outfile)


if __name__ == "__main__":
	main(sys.argv)
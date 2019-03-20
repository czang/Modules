#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

def cutrawdata(infile, number, outfile):
	file = open(infile,'r')
	output = open(outfile,'w')
	i = 0
	for line in file:
		if i < number:
			output.write(line)
		i+=1
	file.close()
	output.close()

def cutrawdata_tail(infile, number, outfile):
	file = open(infile,'r')
	output = open(outfile,'w')
	i = 0
	for line in file:
		if i >= number:
			output.write(line)
		i+=1
	file.close()
	output.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-n", "--threshold", action="store", type="int", dest="number", help="tag number in outfile", metavar="<int>")
	parser.add_option("-t", "--type", action="store", type="int", dest="worktype", help="1 for from head, -1 for from tail", metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	if opt.worktype ==1:
		cutrawdata(opt.infile, opt.number, opt.outfile)
	elif opt.worktype ==-1:
		cutrawdata_tail(opt.infile, opt.number, opt.outfile)
	else:
		print "Type error..."


if __name__ == "__main__":
	main(sys.argv)
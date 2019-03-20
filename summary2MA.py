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
"""


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="summary", metavar="<file>", help="summary file")
	parser.add_option("-a", "--a", action="store", type="float", dest="factora", metavar="<float>", help="total tag count for A")
	parser.add_option("-b", "--b", action="store", type="float", dest="factorb", metavar="<float>", help="total tag count for B")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	o = open(opt.output,'w')
	f = open(opt.summary, 'r')
	i = 1
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			avga = log(max(atof(sline[3]),0.5)/opt.factora * 1000000,2)
			avgb = log(max(atof(sline[4]),0.5)/opt.factorb * 1000000,2)
			avg = avga + avgb
			maxa = max(avga, avgb) * 2
			o.write(sline[0] + '\t' + sline[1] + '\t' + sline[2] + '\t' + str(i) + '\t' + sline[5] + '\t' + str(avg) + '\t' + str(maxa) + '\n')
			i += 1
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
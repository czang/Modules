#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--genefile", action="store", type="string", dest="infile", metavar="<file>", help="chromosome name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.infile+'.txt', 'r')
	genes = f.readline()
	genes = genes.strip()
	f.close()
	o = open('nearestpeak_'+opt.infile+'.sh', 'w')
	o.write('#!/bin/sh'+'\n')
	o.write('for GENE in ' + genes + '\n')
	o.write('do\n')
	o.write('\tsh nearest_peak.sh ' + '$GENE\n')
	o.write('done\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)
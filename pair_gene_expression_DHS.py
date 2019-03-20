#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *
import utility

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile1", action="store", type="string", dest="inputfile", metavar="<file>", help="file name except gene id prefixes")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="listfile", metavar="<file>", help="file of complete gene list with expression")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	g = open(opt.listfile, 'r')
	o = open(opt.out_file, 'w')
	i = 0
	for line in g:
		if (not re.match("track", line)) and (not re.match("variableStep", line)) and (not re.match("#", line)):
			line = line.strip()
			sline = line.split()
			gene = sline[0]
			f = sline[0] + opt.inputfile
			if utility.fileExists(f):
				i += 1
				print i
				p = open(f, 'r')
				for line in p:
					pline = line.strip()
					o.write(gene + '\t' + pline + '\t' + '\t'.join(sline[1:]) + '\n')
				p.close()
	o.close()
	g.close()
	

if __name__ == "__main__":
	main(sys.argv)

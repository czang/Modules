#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import gene_set_manipulation


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();




def main(argv):
	parser = OptionParser()
	parser.add_option("-g", "--gene", action="store", type="string", dest="id", metavar="<file>", help="gene ID")
	parser.add_option("-m", "--matrix", action="store", type="string", dest="matrix_file", metavar="<file>", help="output DNase matrix prefix")
	parser.add_option("-c", "--coefficient", action="store", type="string", dest="coeff", metavar="<file>", help="mlr coefficient file prefix")
	parser.add_option("-d", "--DNase", action="store", type="string", dest="dnasefile", metavar="<file>", help="DNase measurement file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	
	peaklist = gene_set_manipulation.get_gene_list(opt.matrix_file+opt.id, 0)
	gene_set_manipulation.output_UCSCsubset_in_file(opt.dnasefile, peaklist, "temp_"+opt.id)
	
	measurelist = []
	measurelist.append(1.0)
	f = open("temp_"+opt.id, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			measurelist.append(atof(sline[5]))
	f.close()
	
	coefficientlist = []
	f = open(opt.coeff+opt.id, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > 1:
				coefficientlist.append(atof(sline[1]))
	f.close()
	
	assert len(measurelist) == len(coefficientlist)
	
	expression = 0.0
	for i in range(0, len(coefficientlist)):
		expression += coefficientlist[i] * measurelist[i]
	
	print opt.id, expression
		

if __name__ == "__main__":
	main(sys.argv)
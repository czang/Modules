#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
#import BED;
#import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


modifications = ['H2AZ_promoter','H33_promoter','double_promoter','input_promoter','H2AZ_genebody','H33_genebody','double_genebody','input_genebody','H2AZ_end','H33_end','double_end','input_end']



def gene_modification(datafile):
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result[sline[0]] = atof(sline[1])
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfiles", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	complete_data = {}
	for mod in modifications:
		complete_data[mod] = gene_modification(mod+opt.infile)
	
	f = open(opt.output,'w')
	f.write('#gene')
	for mod in modifications:
		f.write('\t' + mod)
	for gene in complete_data['H2AZ_promoter'].keys():
		f.write('\n'+gene)
		for mod in modifications:
			if complete_data[mod][gene] > 0.1:
				f.write('\t1')
			else:
				f.write('\t0')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
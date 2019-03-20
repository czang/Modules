#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();



def make_data_dic(infile):
	file = open(infile,'r')
	result = {}
	x = {}
	y = {}
	for line in file:
		line = line.strip()
		sline = line.split()
		x[sline[0]] = atof(sline[1])
		y[sline[0]] = atof(sline[2])
	result['x'] = x
	result['y'] = y
	return result

def data_subtraction(a,b):
	result = {}
	x = {}
	y = {}
	for gene in a['x'].keys():
		assert a['y'].has_key(gene) and b['x'].has_key(gene) and b['y'].has_key(gene) and a['x'][gene] == b['x'][gene]
		x[gene] = a['x'][gene]
		y[gene] = a['y'][gene] - b['y'][gene]
	result['x'] = x
	result['y'] = y
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--datafile1", action="store", type="string", dest="datafile1", metavar="<file>", help="file with H3K4me3 correlation data")
	parser.add_option("-b", "--datafile2", action="store", type="string", dest="datafile2", metavar="<file>", help="file with H3K27me3 correlation data")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output data file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	result = data_subtraction(make_data_dic(opt.datafile1),make_data_dic(opt.datafile2))
	f = open(opt.output,'w')
	for gene in result['x'].keys():
		f.write(gene + '\t' + str(result['x'][gene]) + '\t' + str(result['y'][gene]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
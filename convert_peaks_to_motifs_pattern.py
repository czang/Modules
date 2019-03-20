#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *

## get BED module
#import BED;
#import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def file_dir(datafile): 
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > 1:
				result[sline[0]] = sline[1:]
	file.close()
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--genepeakfile", action="store", type="string", dest="genefile", metavar="<file>", help="input gene peak file name")
	parser.add_option("-m", "--peakmotiffile", action="store", type="string", dest="peakmotif", metavar="<file>", help="peak motif file name")
	parser.add_option("-l", "--patternlist", action="store", type="string", dest="patterns", metavar="<file>", help="list of all pattern variables")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	patternlist = get_gene_list(opt.patterns, 0)
	motif_dir = file_dir(opt.peakmotif)
	
	infile = open(opt.genefile, 'r')
	outfile = open(opt.output, 'w')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > 1:
				List = []
				for peak in sline[1:]:
					if peak in motif_dir.keys():
						candidates = motif_dir[peak]
						if len(candidates) > 0:
							for item in candidates:
								if not item in List:
									List.append(item)
				if len(List) > 0:
					pattern = []
					for item in patternlist:
						if item in List:
							pattern.append('1')
						else:
							pattern.append('0')
					outfile.write(sline[0] + '\t' + '\t'.join(pattern) + '\n')
				else:
					outfile.write(sline[0] + '\t' + '\t'.join(['0']*len(patternlist)) + '\n')
	infile.close()
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)
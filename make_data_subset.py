#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
from gene_set_manipulation import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--datafile", action="store", type="string", dest="datafile", metavar="<file>", help="data file with gene info in the first column")
	parser.add_option("-b", "--genelist", action="store", type="string", dest="genelist", metavar="<file>", help="file with gene list in the first column")
	parser.add_option("-c", "--directory1", action="store", type="string", dest="indir", metavar="<file>", help="input file directory")
	parser.add_option("-d", "--directory2", action="store", type="string", dest="outdir", metavar="<file>", help="output file directory")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	gene_list = get_gene_list(opt.genelist,0)
	output_UCSCsubset_in_file (opt.indir+'/'+opt.datafile, gene_list, opt.outdir+'/'+opt.datafile)
	print 'Result written in', str(opt.outdir+'/'+opt.datafile)



if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
from species_chroms import *
from associate_island_with_genes import *

Dir = os.getcwd();

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="islands file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-f", "--outfilemarkedislands", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for occupied genes")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
		
	marked_islands = find_genic_islands(opt.bedfile, opt.species, opt.known_genes, opt.promoter_extension, opt.region_type)
	output_bed_graph(marked_islands, opt.out_file, opt.species);	

	
if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
from gene_set_manipulation import *
from associate_island_with_genes import *
from species_chroms import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--infile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="input bed summary graph file")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="input bed summary graph file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0);
	bed_vals_2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0);
	for chrom in bed_vals_1.keys():
		if chrom not in chroms:
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	for chrom in bed_vals_2.keys():
		if chrom not in chroms:
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	shared_number = 0
	onlyone_number = 0
	onlytwo_number = 0
	for chrom in bed_vals_1.keys():
		if chrom in bed_vals_2.keys():
			bed_graph_list_1 = bed_vals_1[chrom]
			bed_graph_list_2 = bed_vals_2[chrom]
			start_list_1 = []
			start_list_2 = []
			for item in bed_graph_list_1:
				start_list_1.append(item.start)
			for item in bed_graph_list_2:
				start_list_2.append(item.start)
			compare = gene_comparison(start_list_1, start_list_2)
			shared_number += len(compare['shared'])
			onlyone_number += len(compare['only in 1'])
			onlytwo_number += len(compare['only in 2'])
		else:
			onlyone_number += len(bed_vals_1[chrom])
	print 'shared:', shared_number, ', file1:', onlyone_number, ', file2:', onlytwo_number


if __name__ == "__main__":
	main(sys.argv)
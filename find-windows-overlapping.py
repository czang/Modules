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
import separate_by_chrom


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--infile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="input bed summary graph file")
	parser.add_option("-b", "--infile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="input bed summary graph file")
	parser.add_option("-p", "--result1", action="store", type="string", dest="result1", metavar="<file>", help="result file name for bed file 1")
	parser.add_option("-q", "--result2", action="store", type="string", dest="result2", help="result file name for bed file 2", metavar="<file>")

	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separate_by_chrom.separateByChrom(chroms, opt.bedfile1, '.bed1')
	separate_by_chrom.separateByChrom(chroms, opt.bedfile2, '.bed2')
	
	for chrom in chroms:
		if separate_by_chrom.fileExists(chrom+'.bed1') and separate_by_chrom.fileExists(chrom+'.bed2'):
			bed_vals_1 = BED.BED(opt.species, chrom+'.bed1', "BED_GRAPH", 0)
			bed_vals_2 = BED.BED(opt.species, chrom+'.bed2', "BED_GRAPH", 0)
			bed_graph_list_1 = bed_vals_1[chrom]
			bed_graph_list_2 = bed_vals_2[chrom]
			start_list_2 = []
			for item in bed_graph_list_2:
				start_list_2.append(item.start)
			overlapped_list_1 = []
			overlapped_list_2 = []
			non_list_1 = []
			non_list_2 = []
			overlapped_starts = []
			for i in range(0, len(bed_graph_list_1)):
				if bed_graph_list_1[i].start in start_list_2:
					overlapped_list_1.append(bed_graph_list_1[i])
					overlapped_starts.append(bed_graph_list_1[i].start)
					
				else:
					non_list_1.append(bed_graph_list_1[i])
			for i in range(0, len(bed_graph_list_2)):
				if bed_graph_list_2[i].start in overlapped_starts:
					overlapped_list_2.append(bed_graph_list_2[i])
				else:
					non_list_2.append(bed_graph_list_2[i])
			file = open(.overlap1)
		else: 



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
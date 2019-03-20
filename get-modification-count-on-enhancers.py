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
import get_total_tag_counts
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

'''This module is to calculate the modification data on any region given as UCSC known gene format. The region can have different-length upstream and downstream extensions. Given a p-value, all the windows having tags more than the window threshold is counted.'''


def get_threshold_table(thresholdfile):
	result_table = {}
	file = open(thresholdfile,'r')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result_table[atoi(sline[0])] = atoi(sline[1]);
	file.close()
	return result_table


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-t", "--startextension", action="store", type="int", dest="start_extension", help="upstream threshold of gene region", metavar="<int>")
	parser.add_option("-e", "--endextension", action="store", type="int", dest="end_extension", help="downstream threshold of gene region", metavar="<int>")
	parser.add_option("-p", "--thresholdfile", action="store", type="string", dest="thresholdfile", help="file with the list of thresholds for different number of windows", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0);
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	threshold_table = get_threshold_table(opt.thresholdfile)
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.place_domains_on_regions_with_changing_thresholds(bed_graph_list, genes, opt.start_extension, opt.end_extension, 200, threshold_table)
			for gene in current_genes.keys():
				occupied_genes[gene] = current_genes[gene]
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
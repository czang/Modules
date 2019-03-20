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

'''This module is to calculate the modification data specially on gene promoters. Modification data are gained from summary graph file. Given a p-value, if total tag count on the promoter region is more than the threshold, the number is recorded with the gene. Otherwise record 0.'''

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
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
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile)
	average = total_tag_counts / genomelength * 2 * float(opt.promoter_extension)
	threshold = associate_island_with_genes.find_threshold(opt.pvalue, average)
	print "The tag threshold on promoter is", threshold;
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			#total += len(genes);
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.place_domains_on_genes(bed_graph_list, genes, 'Promoter', opt.promoter_extension);
			for gene in current_genes.keys():
				if current_genes[gene] > threshold:
					total +=1
	print opt.bedfile, total_tag_counts, float(total)/12541.0


if __name__ == "__main__":
	main(sys.argv)
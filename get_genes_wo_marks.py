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
import gene_set_manipulation
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();




def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="islands file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-f", "--outfileoccupiedgenes", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for occupied genes")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter'];
	if opt.region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);
		
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0); #islands 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;

	for chrom in chroms:
		genes = coords[chrom];
		total += len(genes);
		bed_graph_list = bed_vals[chrom];
		current_genes = associate_island_with_genes.place_domains_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension);
		for gene in current_genes.keys():
			if current_genes[gene] <.001:
				occupied_genes[gene] = current_genes[gene];
	print "Number of marked genes is: ", len(occupied_genes), " out of total of ", total, " genes in "+ opt.species ;
	
	gene_list = occupied_genes.keys();
	gene_set_manipulation.output_UCSCsubset_in_file(opt.known_genes, gene_list, opt.out_file);
	
if __name__ == "__main__":
	main(sys.argv)

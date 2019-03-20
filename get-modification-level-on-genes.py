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


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="windowsize", help="window size of summary graph file, in bps", metavar="<int>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
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
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0);
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile)
	average = total_tag_counts / genomelength * float(opt.windowsize)
	threshold = associate_island_with_genes.find_threshold(opt.pvalue, average) # calculating threshold
	pseudocount = threshold / 2 + 1;
	print "The window threshold of tags is", threshold;
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.place_domain_windows_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension, threshold);
			for gene in current_genes.keys():
				if current_genes[gene] > .1:
					occupied_genes[gene] = float(current_genes[gene])/float(total_tag_counts);
				else:
					occupied_genes[gene] = float(pseudocount)/float(total_tag_counts)
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
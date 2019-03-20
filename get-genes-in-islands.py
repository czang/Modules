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
import gene_set_manipulation

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island bed file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
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
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED3", 0); 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	inside_number = 0
	overlap_number = 0
	
	distance_list = []
	
	inside_list = []
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			bed_graph_list = bed_vals[chrom];
			(inside, overlap, other) = associate_island_with_genes.identify_gene_on_islands(bed_graph_list, genes, opt.region_type, opt.promoter_extension);
			inside_list += inside
			inside_number += len(inside)
			overlap_number += len(overlap)
			current_genes = associate_island_with_genes.find_closest_island_distance_around_gene(bed_graph_list, other);
			for gene in current_genes.keys():
				occupied_genes[gene] = current_genes[gene]
				distance_list.append(atoi(current_genes[gene]))
	
	distance_list.sort()
	N10K = bisect.bisect_left(distance_list, 10000)
	N20K = bisect.bisect_left(distance_list, 20000) - bisect.bisect_left(distance_list, 10000)
	N30K = bisect.bisect_left(distance_list, 30000) - bisect.bisect_left(distance_list, 20000)
	N40K = bisect.bisect_left(distance_list, 40000) - bisect.bisect_left(distance_list, 30000)
	N50K = bisect.bisect_left(distance_list, 50000) - bisect.bisect_left(distance_list, 40000)
	Nfar = len(distance_list) - bisect.bisect_left(distance_list, 50000)
	
	print inside_number, overlap_number, N10K, N20K, N30K, N40K, N50K, Nfar
	gene_set_manipulation.output_UCSCsubset_in_file (opt.known_genes, inside_list, opt.out_file)
	#f = open(opt.out_file, 'w')
	#for g in occupied_genes.keys():
	#	f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	#f.close()


if __name__ == "__main__":
	main(sys.argv)
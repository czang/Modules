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

'''This module is to calculate the modification data specially on gene promoters. Modification data are gained from summary graph file. Given a p-value, if total tag count on the promoter region is more than the threshold, the number is recorded with the gene. Otherwise record 0. 
The threshold calculation is taken into account a background library, known as IgG or input library. The input background file is a list of genes with the significant counts on the SAME promoter region (i.e. promoter extension must be same), the genelist must have the same size of the all genes to be caculated. The average tag counts on a promoter is the maximum value between the real average distributed value and the rescaled background distributed value.'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-b", "--background_file", action="store", type="string", dest="background", metavar="<file>", help="file with control scores for all genes")
	parser.add_option("-i", "--backgroundcounts", action="store", type="int", dest="backgroundtotal", help="total tag counts of background library", metavar="<int>")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
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
	print "Average tag counts on promoter is", average
	
	input_count_dic = gene_set_manipulation.get_gene_float_dic(opt.background, 1)
	threshold_dic = {}
	for gene in input_count_dic.keys():
		rescaled_average = max(input_count_dic[gene] * total_tag_counts / opt.backgroundtotal, average)
		threshold_dic[gene] = associate_island_with_genes.find_threshold(opt.pvalue, rescaled_average)
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			total += len(genes);
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.place_domains_on_genes(bed_graph_list, genes, 'Promoter', opt.promoter_extension);
			for gene in current_genes.keys():
				if gene in threshold_dic.keys():
					if current_genes[gene] > threshold_dic[gene]:
						occupied_genes[gene] = int(current_genes[gene])
					else:
						occupied_genes[gene] = 0
				else: 
					print gene, "does not have a corresponded threshold."
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
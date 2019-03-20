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


def result_comparison(Result_1, Result_2, pseudocount, normalization):
	result = {}
	if len(Result_1.keys()) != len(Result_2.keys()):
		print 'result error'
	for name in Result_1.keys():
		if Result_2.has_key(name):
			if Result_1[name] < 0.01:
				Result_1[name] += pseudocount
				Result_2[name] += pseudocount
			else:
				Result_1[name] = float(Result_1[name]) / float(normalization)
				if Result_2[name] < 0.01:
					Result_1[name] += pseudocount
					Result_2[name] += pseudocount
			result[name] = log(float(Result_1[name]) / float(Result_2[name]),10)
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="summary graph file for intial state")
	parser.add_option("-b", "--summarygraphfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="summary graph file for final state")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="windowsize", help="window size of summary graph file, in bps", metavar="<int>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 18:
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
	
	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0);
	bed_vals_2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0);
	#for chrom in chroms:
		#if (chrom not in bed_vals_1.keys()) or (chrom not in bed_vals_2.keys()):
			#print chrom, " name is not the same as the stored one";
			#sys.exit(1);
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	total_counts_1 = float(get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile1))
	total_counts_2 = float(get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile2))
	threshold1 = associate_island_with_genes.find_threshold(opt.pvalue, total_counts_1 / genomelength * float(opt.windowsize))
	threshold2 = associate_island_with_genes.find_threshold(opt.pvalue, total_counts_2 / genomelength * float(opt.windowsize))
	pseudocount = (threshold1 + threshold2) / 4 + 1;
	print "The window thresholds of tags are", threshold1, 'and', threshold2;
	
	coords = UCSC.KnownGenes(opt.known_genes);
	#occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;
	Result = {}
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			total += len(genes);
			bed_graph_list_1 = bed_vals_1[chrom];
			current_genes_1 = associate_island_with_genes.place_domain_windows_on_genes(bed_graph_list_1, genes, opt.region_type, opt.promoter_extension, threshold1);
			bed_graph_list_2 = bed_vals_2[chrom];
			current_genes_2 = associate_island_with_genes.place_domain_windows_on_genes(bed_graph_list_2, genes, opt.region_type, opt.promoter_extension, threshold2);
			Result[chrom] = result_comparison(current_genes_2, current_genes_1, pseudocount, total_counts_2 / total_counts_1)
	
	f = open(opt.out_file, 'w')
	for chrom in Result.keys():
		for g in Result[chrom].keys():
			f.write(g + '\t' + str(Result[chrom][g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)

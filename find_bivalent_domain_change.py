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


def result_comparison(Result_1, Result_2, total_1, total_2, pseudocount):
	result = {}
	if len(Result_1.keys()) != len(Result_2.keys()):
		print 'result error'
	for name in Result_1.keys():
		if Result_2.has_key(name):
			if Result_1[name] < 0.01 and Result_2[name] < 0.01:
				result[name] = 1.0
			elif Result_1[name] < 0.01:
				result[name] = Result_2[name] / total_2 / pseudocount * total_1
			elif Result_2[name] < 0.01:
				result[name] = Result_1[name] / total_1 / pseudocount * total_2
			else:
				result[name] = Result_1[name] / total_1 / Result_2[name] * total_2
			if result[name] < 1.0:
				result[name] = 1.0 / result[name]
	return result


def result_combination(Result_1, Result_2):
	result = {}
	if len(Result_1.keys()) != len(Result_2.keys()):
		print 'result error'
	for g in Result_1.keys():
		if Result_2.has_key(g):
			result[g] = Result_1[g] * Result_2[g]
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="K4file1", metavar="<file>", help="H3K4me summary of celltype 1")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="K27file1", metavar="<file>", help="H3K27me summary of celltype 1")
	parser.add_option("-c", "--bedfile3", action="store", type="string", dest="K4file2", metavar="<file>", help="H3K4me summary of celltype 2")
	parser.add_option("-d", "--bedfile4", action="store", type="string", dest="K27file2", metavar="<file>", help="H3K27me summary of celltype 2")
	parser.add_option("-k", "--bivalent_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known bivalent genes in UCSC format")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
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
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	
	K4file1_total = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.K4file1)
	K4file2_total = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.K4file2)
	K27file1_total = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.K27file1)
	K27file2_total = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.K27file2)
	
	K4_bed_vals_1 = BED.BED(opt.species, opt.K4file1, "BED_GRAPH", 0);
	K4_bed_vals_2 = BED.BED(opt.species, opt.K4file2, "BED_GRAPH", 0);
	K27_bed_vals_1 = BED.BED(opt.species, opt.K27file1, "BED_GRAPH", 0);
	K27_bed_vals_2 = BED.BED(opt.species, opt.K27file2, "BED_GRAPH", 0);
	
	K4_threshold1 = associate_island_with_genes.find_threshold(0.001, K4file1_total / genomelength * 200.0)
	K4_threshold2 = associate_island_with_genes.find_threshold(0.001, K4file2_total / genomelength * 200.0)
	K27_threshold1 = associate_island_with_genes.find_threshold(0.001, K27file1_total / genomelength * 200.0)
	K27_threshold2 = associate_island_with_genes.find_threshold(0.001, K27file2_total / genomelength * 200.0)
	pseudocount = (K4_threshold1 + K4_threshold2 + K27_threshold1 + K27_threshold2) / 8 + 1;
	print "The window thresholds of tags:", K4_threshold1, K27_threshold1, K4_threshold2, 'and', K27_threshold2;
	
	coords = UCSC.KnownGenes(opt.known_genes);
	#occupied_genes = {}; # a BED object of BED_GRAPH elements
	Result = {}
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			K4_list_1 = K4_bed_vals_1[chrom];
			K4_genes_1 = associate_island_with_genes.place_domain_windows_on_genes(K4_list_1, genes, opt.region_type, opt.promoter_extension, K4_threshold1);
			K4_list_2 = K4_bed_vals_2[chrom];
			K4_genes_2 = associate_island_with_genes.place_domain_windows_on_genes(K4_list_2, genes, opt.region_type, opt.promoter_extension, K4_threshold2);
			K27_list_1 = K27_bed_vals_1[chrom];
			K27_genes_1 = associate_island_with_genes.place_domain_windows_on_genes(K27_list_1, genes, opt.region_type, opt.promoter_extension, K27_threshold1);
			K27_list_2 = K27_bed_vals_2[chrom];
			K27_genes_2 = associate_island_with_genes.place_domain_windows_on_genes(K27_list_2, genes, opt.region_type, opt.promoter_extension, K27_threshold2);
			Result[chrom] = result_combination(result_comparison(K4_genes_1, K4_genes_2, K4file1_total, K4file2_total, pseudocount), result_comparison(K27_genes_1, K27_genes_2, K27file1_total, K27file2_total, pseudocount))
	f = open(opt.out_file,'w')
	for chrom in Result.keys():
		for g in Result[chrom].keys():
			f.write(g + '\t' + str(Result[chrom][g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)

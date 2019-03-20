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
#import get_total_tag_counts
import SeparateByChrom
import associate_island_with_genes
import associate_tags_with_regions
import Utility

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def find_index_list_w_max(name_list, count_list):
	assert len(name_list) == len(count_list)
	out_list = []
	i = 0
	j = i + 1
	while i < len(name_list):
		current = i
		j = i + 1
		while j < len(name_list):
			if name_list[j] == name_list[i]:
				if  count_list[j] > count_list[current]:
					current = j
				j += 1
			else:
				break
		out_list.append(current)
		i = j
	return out_list


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-f", "--fragmentsize", action="store", type="int", dest="fragment_size", help="fragment size of ChIP-seq reads, in bps", metavar="<int>")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--'Promoter', 'GeneBody', 'GenePromoter', or 'GeneExtended'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>")
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
	
	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter', 'GeneExtended'];
	if opt.region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);
	
	if Utility.fileExists(opt.bedfile):
		SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1');
	else:
		print opt.bedfile, " not found";
		sys.exit(1)
	
	#genomelength = 0.0
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	#total = 0;
	
	outputfile = open(opt.out_file, 'w')
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			(gene_name_list, region_start_list, region_end_list) = associate_island_with_genes.get_gene_feature_region_lists(coords[chrom], opt.region_type, opt.promoter_extension)
			
			tag_position_list = []
			read_file = chrom + ".bed1";
			f = open(read_file,'r')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_position_list.append(associate_tags_with_regions.tag_position(sline, opt.fragment_size))
			f.close();
			tag_position_list.sort()
			tag_count_list = associate_tags_with_regions.find_readcount_on_regions(tag_position_list, region_start_list, region_end_list)
			
			assert len(gene_name_list) == len(tag_count_list)
			index_list = find_index_list_w_max(gene_name_list, tag_count_list)
			
			for i in index_list:
				g = genes[i]
				outputfile.write(g.name + "\t" + chrom + "\t" + g.strand + "\t" + str(g.txStart) + "\t" + str(g.txEnd) + "\t" + str(tag_count_list[i]) + "\t0\t0\t0\t0\n")
	outputfile.close()
	SeparateByChrom.cleanup(chroms, '.bed1');


if __name__ == "__main__":
	main(sys.argv)
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

'''This module is to calculate the modification data specially on gene 3' ends with different extended lengths. Modification data are gained from summary graph file. Given a p-value, if total tag count on the promoter region is more than the threshold, the number is recorded with the gene. Otherwise record 0.'''

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-t", "--startextension", action="store", type="int", dest="start_extension", help="upstream threshold of gene region", metavar="<int>")
	parser.add_option("-e", "--endextension", action="store", type="int", dest="end_extension", help="downstream threshold of gene region", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
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
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile)
	average = total_tag_counts / genomelength * float(opt.start_extension + opt.end_extension)
	#threshold = associate_island_with_genes.find_threshold(opt.pvalue, average)
	#print "The tag threshold on promoter for", opt.bedfile, ':', threshold;
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;
	
	Both_count = 0
	Single_count = 0
	No_count = 0
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			total += len(genes);
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.place_domains_on_promoters(bed_graph_list, genes, opt.start_extension, opt.end_extension);
			end_genes = associate_island_with_genes.place_domains_on_gene_ends(bed_graph_list, genes, opt.end_extension, opt.start_extension);
			assert len(current_genes) == len(end_genes)
			for gene in current_genes.keys():
				if gene in end_genes.keys():
					if current_genes[gene] > 0.01 and end_genes[gene] > 0.01:
						occupied_genes[gene] = 2
						Both_count += 1
					elif current_genes[gene] > 0.01:
						occupied_genes[gene] = 1
						Single_count += 1
					elif end_genes[gene] > 0.01:
						occupied_genes[gene] = 1
						Single_count += 1
					else:
						occupied_genes[gene] = 0
						No_count += 1
	
	print "Islands found in both TSS and TES: ", Both_count
	print "Islands found in one side bound:", Single_count
	print "No islands found in boudaries:", No_count
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
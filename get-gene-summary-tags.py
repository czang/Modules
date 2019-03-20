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

def place_domains_on_genes(bed_graph_list, genes, region_type, promoter_extension, threshold):
	"""
	region_type:
	*Promoter: TSS promoter region
	*GeneBody:  GB gene body sbtract the part already in TSS,  mutually exclusive with TSS
	*GenePromoter: G GB+TSS
	
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of gene name: tag_number
	"""
	occupied_gene={};
	if len(bed_graph_list) != 0:
		island_start_list=[];
		island_end_list=[];
		island_value_list =[];
		for item in bed_graph_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
			island_value_list.append(item.value);
		for g in genes:
			# gene body and TSS are mutually exclusive
			if plus.match(g.strand):
				TSS_start = g.txStart - promoter_extension
				TSS_end = g.txStart + promoter_extension
				GB_start = TSS_end
				GB_end = g.txEnd
				G_start = TSS_start
				G_end = GB_end
			# When strand is negative, the regions are still sequential
			elif minus.match(g.strand):
				TSS_start = g.txEnd - promoter_extension
				TSS_end = g.txEnd + promoter_extension
				GB_start = g.txStart
				GB_end = TSS_start
				G_start = GB_start
				G_end = TSS_end
			else: print "Gene strand identification crisis!";
			
			if region_type == 'Promoter':
				start = TSS_start
				end = TSS_end
			elif region_type == 'GeneBody': # start could be bigger than end if TSS is chosen big.
				start = GB_start
				end = GB_end
			elif region_type == 'GenePromoter':
				start = G_start
				end = G_end
			tag = 0;
			association_list = associate_island_with_genes.place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			#if g.name == '207480_s_at':
				 ##print start, end, association_list;
				 ##print island_start_list[18], island_end_list[18];
				 #print island_start_list[19], island_end_list[19];
				 #print island_start_list[20], island_end_list[20];
				 #print island_start_list[21], island_end_list[21];
			for i in association_list: 
				if island_value_list[i] >= threshold:
					tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0;
	return occupied_gene;


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="islands file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-t", "--windowthreshold", action="store", type="float", dest="windowthreshold", help="window shreshold, window that have no less than this number of tags are counted", metavar="<int>")
	parser.add_option("-b", "--backgroundlevel", action="store", type="float", dest="zerocount", help="background score for genes that there is no island on it", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
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
		current_genes = place_domains_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension, opt.windowthreshold);
		for gene in current_genes.keys():
			if current_genes[gene] >.01:
				occupied_genes[gene] = current_genes[gene];
			else:
				occupied_genes[gene] = opt.zerocount
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()
	
	#gene_list = occupied_genes.keys();
	#gene_set_manipulation.output_UCSCsubset_in_file(opt.known_genes, gene_list, opt.out_file);
	
if __name__ == "__main__":
	main(sys.argv)

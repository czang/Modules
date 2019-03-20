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


def get_total_tags_on_islands(islands):
	'''return the average tag density of the islands in one chromosome. It is the ratio between the total tag number in all islands and the total island lengths '''
	total_counts = 0.0
	for island_index in range(0, len(islands)):
		total_counts = total_counts + islands[island_index].value;
	return float(total_counts);

def get_total_tags_on_islands_genome(bed_vals):
	total = 0.0;
	for chrom in bed_vals.keys():
		islands =  bed_vals[chrom];	
		total += get_total_tags_on_islands(islands);
	return total;

def get_normalization_factor(bed_val_1, bed_val_2):
	'''
	Because different modification has quite different number of tags, the absolution tag density needs to be
	adjusted. It is not clear how best to define this as the two different modification might have totally different
	characteristics.
	'''
	total_1 = get_total_tags_on_islands_genome(bed_val_1);
	total_2 = get_total_tags_on_islands_genome(bed_val_2);
	#print total_1, total_2, total_1/total_2;
	return total_1/total_2;
	
def find_bivalent_islands(bed_vals_1, bed_vals_2, extention, chrom, bed_val_2_normalization_factor):	
	'''
	Return the list of bivalent domains, each element in the list is a BED_GRAPH object.
	Use K4 for bed_val_1

	For each island i in bed_vals_1, find if there are any islands in bed_vals_2 that is within extension
	of it. If yes, this island i is in a bivalent domain and will be out put in a fil f. The information
	writen to f is
	the start of the bivalent domain
	the end of the bivalent domain
	the effective tag counts on the bivalent domain

	The bivalent domain is defined in terms of bed_val_1, it is the union of the island i and those islands in
	bed_val_2 that are associated with it. What might be better is to find the cluster of islands for both modification1 and modification2. Not sure how to make it happen yet.

	What is needed is the pairwise information for use of comparison with other cell-types. We are using K4 islands as the coordinate for the bivalent domain.
	'''
	bivalent_domains =[];
	
	if chrom in bed_vals_1.keys() and chrom in bed_vals_2.keys():
		islands_1 = bed_vals_1[chrom]
		islands_2 = bed_vals_2[chrom]
		
		if len(islands_1) != 0 and len(islands_2) != 0:
			islands_2_start = bed_vals_2.getStarts(chrom);
			islands_2_end = bed_vals_2.getEnds(chrom);
				
			for island_index in range(0, len(islands_1)):# run through each island in islands_1
				"""
				classic binary search --
				for bisect_left: If 'search' is in 'list', its index is returned.  If not, the index where 'search' should be in 'list' is returned
				for bisect_right:If 'search' is in 'list', its index+1 is returned. If not, the index where 'search' should be in 'list' is returned

				Think of bisect as looking for the position to add the 'search' into the already sorted list. 

				index_start = bisect.bisect_left(islands_2_end, islands_1[island_index].start - extention):
				Search islands_1's start in list islands_2_end
				If islands_1[island_index] overlaps with islands_2, the first one would be islands_2[index_start]

				index_end = bisect.bisect_right(islands_2_start, islands_1[island_index].end + extention)
				Search islands_1's end in list islands_2_start.
				If islands_1[island_index] overlaps with islands_2, the last one would have index index_end-1. 
				"""
				index_start = bisect.bisect_left(islands_2_end, islands_1[island_index].start - extention)
				index_end = bisect.bisect_right(islands_2_start, islands_1[island_index].end + extention)
					
				if index_end > index_start:  # Overlap found
					bivalent_region_start = islands_1[island_index].start;
					bivalent_region_end = islands_1[island_index].end;
					
					#print island_index, index_start, index_end;
					#bivalent_region_start = min(islands_1[island_index].start , islands_2[index_start].start);
					#bivalent_region_end = max(islands_1[island_index].end, islands_2[index_end-1].end);
					tag_counts_in_islands_2 = 0;
					# Find all the islands in islands_2 that overlaps with islands_1[island_index]
					for i in range(index_start, min(index_end, len(islands_2_end))):
						tag_counts_in_islands_2 += islands_2[i].value;
					bivalent_region_value = int(islands_1[island_index].value + bed_val_2_normalization_factor * tag_counts_in_islands_2);
					bed = BED.BED_GRAPH(chrom, bivalent_region_start, bivalent_region_end, bivalent_region_value);
					bivalent_domains.append(bed);

					#f.write(chrom + ' ' + str(bivalent_region_start) + ' ' + str(bivalent_region_end) + ' ' + str(bivalent_region_value)+'\n');
			#bivalent_domains.sort(key=operator.attrgetter('value'));
			print chrom, 'total island number   1:',len(islands_1),', 2:',len(islands_2),'. Bivalent domains with respect to 1:', len(bivalent_domains), str(int(float(100*len(bivalent_domains))/float(len(islands_1))))+'%';
		else:
			print chrom, 'total island number   1:',len(islands_1),', 2:',len(islands_2),'. No bivalent domains.'
	else:
		print "The chromosome " + chrom + " is not in this species! ";

	return bivalent_domains;

	

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="islands from modification 1")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="islands from modification 2")
	parser.add_option("-e", "--extention", action="store", type="int", dest="extention",  metavar="<int>", help="distance allowed between bivalent islands")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_gene_file", help="output file name for bivalent genes in UCSC format", metavar="<file>")
	
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

	total_islands_1 = 0;
	bivalent_domains = {};
       	total_domains = 0;

	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0) #islands from modification 1
	bed_vals_2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0) #islands from modification 2

	bed_val_2_normalization_factor = get_normalization_factor(bed_vals_1, bed_vals_2);
	
	coords = UCSC.KnownGenes(opt.known_genes)
	
	for chrom in chroms:
		number_of_islands = len(bed_vals_1[chrom]);
		bivalent_domains[chrom] = find_bivalent_islands(bed_vals_1, bed_vals_2, opt.extention, chrom, bed_val_2_normalization_factor);
		total_islands_1 += len(bed_vals_1[chrom]);
		total_domains += len(bivalent_domains[chrom]);
		#print chrom, number_of_islands, len(bivalent_domain[chrom]);
	print "There are ", total_domains, " bivalent domains out of ", total_islands_1, "islands in " + opt.bedfile1;	
	#associate_island_with_genes.output_bed_graph(bivalent_domains, opt.out_bivalentdomain_file, opt.species);
	
	#file = open(opt.out_gene_file,'w');
	bivalent_gene_list = []
	#overlapped_domains = 0;
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			bed_graph_list = bivalent_domains[chrom];
			occupied_genes = associate_island_with_genes.place_domains_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension);
			for gene in occupied_genes.keys():
				if occupied_genes[gene] > 0.1:
					bivalent_gene_list.append(gene);
	gene_set_manipulation.output_UCSCsubset_in_file(opt.known_genes, bivalent_gene_list, opt.out_gene_file)
		#for k,v in occupied_genes.items(): 
			#if v != 0: file.write(k+" "+str(v)+"\n");
		#overlapped_islands =associate_island_with_genes.find_domains_not_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension);
		#overlapped_domains += len([elem for elem in overlapped_islands if elem!=0]);
	#print overlapped_domains, " bivalent domains are inside the " + opt.region_type + " region";


if __name__ == "__main__":
	main(sys.argv)

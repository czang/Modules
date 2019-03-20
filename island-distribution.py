#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Authors: Chongzhi Zang and Xiaole Shirley Liu
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# czang@jimmy.harvard.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

## get BED module
import BED_MACS;
import UCSC;
from GenomeData import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is used to associate the islands distribution with genes. Of all islands through the genome of a cell, it tells: 
1) How many are in the promotor region 
2) How many are in the gene body region 
3) In gene body region, how many are in exons and how many are in introns, respect with the bp length of exons and introns 
4) How many are in the intergenic region

'''


def find_island_summit(island):
	return int(island.start + island.summit - 1)


def is_location_in_gene(position, genes):
	flag = 0
	index = 0
	outindex = -1
	while flag == 0 and index < len(genes):
		if genes[index].txStart <= position and genes[index].txEnd >= position:
			flag = 1
			outindex = index
		else:
			index += 1
	if outindex == -1:
		return 0
	else:
		return genes[outindex]


def get_island_regions(bed_vals, coords, chrom, TSSthreshold):
	if chrom in bed_vals.keys() and chrom in coords.keys():
		islands = bed_vals[chrom]
		
		genes = coords[chrom]
		TSS_List = [] ##make TSS list
		plus_TSS_list = []
		minus_TSS_list = []
		
		for g in genes: 
			if plus.match(g.strand):
				TSS_List.append(g.txStart);
				plus_TSS_list.append(g.txStart)
			elif minus.match(g.strand):
				TSS_List.append(g.txEnd)
				minus_TSS_list.append(g.txEnd)
		TSS_List.sort()
		plus_TSS_list.sort()
		minus_TSS_list.sort()
		
		Summit = 0
		TSSindex = 0
		
		promoter_islands = []
		exon_islands = []
		intron_islands = []
		intergenic_islands = []
		upstream_distances = []
		
		for island in islands:
			Summit = find_island_summit(island)
			TSSindex = bisect.bisect_left(TSS_List, Summit)
			target_gene = is_location_in_gene(Summit, genes)
			if abs(Summit - TSS_List[min(TSSindex,(len(TSS_List)-1))]) < TSSthreshold or abs(Summit - TSS_List[max(0,(TSSindex-1))]) < TSSthreshold:
				#TSS_Number += 1; 
				promoter_islands.append(island)
			elif target_gene != 0: 
				#GeneBody_Number += 1;
				#print genes[index_end].name
				exon_Start_List = []
				exon_End_List = []
				exon_Starts = target_gene.exonStarts.replace(',',' ')
				exon_Ends = target_gene.exonEnds.replace(',',' ')
				exon_Starts = exon_Starts.split()
				exon_Ends = exon_Ends.split()
				for item in exon_Starts:
					exon_Start_List.append(int(item))
				for item in exon_Ends:
					exon_End_List.append(int(item))
				exon_index_start = bisect.bisect_left(exon_Start_List, Summit)
				exon_index_end = bisect.bisect_left(exon_End_List, Summit)
				if exon_index_start != exon_index_end:
					#exon_Number += 1
					exon_islands.append(island)
				else:
					#intron_Number += 1
					intron_islands.append(island)
			else: 
				#Intergenic_Number += 1 
				intergenic_islands.append(island)
				plus_index = bisect.bisect_left(plus_TSS_list, Summit)
				if plus_index < len(plus_TSS_list):
					plus_distance = abs(Summit - plus_TSS_list[plus_index])
				else:
					plus_distance = 10000000000
				minus_index = bisect.bisect_left(minus_TSS_list, Summit)
				if minus_index > 0:
					minus_distance = abs(Summit - minus_TSS_list[minus_index - 1])
				else:
					minus_distance = 10000000000
				upstream_distances.append(min(plus_distance, minus_distance))
	return promoter_islands, exon_islands, intron_islands, intergenic_islands, upstream_distances




def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="file with islands info in BED MACS format")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-t", "--TSS_threshold", action="store", type="int", dest="TSSthreshold", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	bed_vals = BED_MACS.BED(opt.species, opt.bedfile, "BED_MACS", 0)
	coords = UCSC.KnownGenes(opt.known_genes)
	total_promoter_islands = {}
	total_exon_islands = {}
	total_intron_islands = {}
	total_intergenic_islands = {}
	total_upstream_distances = {}
	promoter_num = 0
	exon_num = 0
	intron_num = 0
	intergenic_num = 0 
	for chrom in chroms: 
		(promoter_islands, exon_islands, intron_islands, intergenic_islands, upstream_distances) = get_island_regions(bed_vals, coords, chrom, opt.TSSthreshold)
		assert len(intergenic_islands) == len(upstream_distances)
		total_promoter_islands[chrom] = promoter_islands
		total_exon_islands[chrom] = exon_islands
		total_intron_islands[chrom] = intron_islands
		total_intergenic_islands[chrom] = intergenic_islands
		total_upstream_distances[chrom] = upstream_distances
		promoter_num += len(promoter_islands)
		exon_num += len(exon_islands)
		intron_num += len(intron_islands)
		intergenic_num += len(intergenic_islands)
	print 'Promoter:', promoter_num, '; Exon:', exon_num, '; Intron:', intron_num, '; Intergenic:',  intergenic_num
	f = open(opt.out_file+'_promoter', 'w')
	for chrom in total_promoter_islands.keys():
		for island in total_promoter_islands[chrom]:
			f.write(chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\t' + str(island.score) + '\n')
	f.close()
	f = open(opt.out_file+'_exon', 'w')
	for chrom in total_exon_islands.keys():
		for island in total_exon_islands[chrom]:
			f.write(chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\t' + str(island.score) + '\n')
	f.close()
	f = open(opt.out_file+'_intron', 'w')
	for chrom in total_intron_islands.keys():
		for island in total_intron_islands[chrom]:
			f.write(chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\t' + str(island.score) + '\n')
	f.close()
	f = open(opt.out_file+'_intergenic', 'w')
	for chrom in total_intergenic_islands.keys():
		for i in range(0, len(total_intergenic_islands[chrom])):
			island = total_intergenic_islands[chrom][i]
			dis = total_upstream_distances[chrom][i]
			f.write(chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\t' + str(island.score) + '\t' + str(dis) + '\n')
	f.close()
	


if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
# Copyright (c) 2009 GWU & NHLBI, NIH
# Authors: Chongzhi Zang, Weiqun Peng and Keji Zhao
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
# zang@gwmail.gwu.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from GenomeData import *
import get_total_tag_counts
import gene_set_manipulation

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def factorial(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;


def factln(m):
	if m<20:
		value = 1.0;
		if m != 0:
			while m != 1:
				value = value*m;
				m = m - 1;
		return log(value);
	else:
		return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;


def poisson(i, average):
	if i<20:
		return exp(-average) * average**i / factorial(i);
	else:
		exponent = -average + i*log(average) - factln(i);
		return exp(exponent);

#def poisson(i, average):
#	return exp(-average) * pow(average,i) / factorial(i);


def find_threshold(pvalue, average):
	value = 1;
	index = 0;
	value -= poisson(index, average);
	while value > pvalue:
		index += 1;
		value -= poisson(index, average);	
	return index;


def get_TSS_list(genes):
	TSS_list = []
	if len(genes) != 0:
		for g in genes:
			if plus.match(g.strand):
				TSS_list.append(g.txStart)
			elif minus.match(g.strand):
				TSS_list.append(g.txEnd)
	return TSS_list


def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def FindIslandCenter(island_start, island_end, StartList, EndList, TagCountList, WindowSize):
	'''Return the average central coordinate of an island, with window counts weighted, i.e. the center of mass of an island. Three lists are from the summary graph'''
	StartPosition = bisect.bisect_left(StartList, island_start)
	EndPosition = bisect.bisect_left(EndList, island_end)
	assert island_start == StartList[StartPosition] and island_end == EndList[EndPosition]
	average = 0.0
	Total = 0.0;
	TotalCounts = 0.0;
	for i in range(StartPosition, EndPosition+1):
		Total = Total + TagCountList[i] * ( StartList[i] + WindowSize / 2 )
		TotalCounts = TotalCounts + TagCountList[i]
	average = Total / TotalCounts
	return int(average)


def generate_island_position_and_value_lists(island_list, summary_list, windowsize):
	'''return two lists, with the centers of mass of islands, and values, respectively'''
	island_position_list = []
	island_value_list = []
	summary_start_list = []
	summary_end_list = []
	summary_value_list = []
	for i in summary_list:
		summary_start_list.append(i.start)
		summary_end_list.append(i.end)
		summary_value_list.append(i.value)
	island_list.sort(key=operator.attrgetter('start'))
	for i in island_list:
		island_position_list.append(FindIslandCenter(i.start, i.end, summary_start_list, summary_end_list, summary_value_list, windowsize))
		island_value_list.append(i.value)
	assert len(island_position_list) == len(island_value_list)
	return (island_position_list, island_value_list)


def place_island_center_on_regions(island_position_list, island_value_list, region_start_list, region_end_list):
	'''returns a list, with total tag number of all islands on this region, order as the region lists'''
	result = []
	assert len(region_start_list) == len(region_end_list)
	assert len(island_position_list) == len(island_value_list)
	for i in range(0, len(region_start_list)):
		tags = 0
		s = bisect.bisect_left(island_position_list, region_start_list[i])
		e = bisect.bisect_right(island_position_list, region_end_list[i])
		if s != e:
			for j in range(s, e):
				tags += island_value_list[j]
		result.append(tags)
	assert len(result) == len(region_end_list)
	return result


def associate_a_location_to_genes(location, start_list):
	"""
	This one compares the location with the starts of the
	genes, and find the closest one for association. 
	
	The algorithm runs through every gene, not particularly efficient

	Return the particular g in genes that associates with this island.
	
	"""
	new_list = []
	for start_position in start_list:
		new_list.append(abs(start_position-location));
	index = new_list.index(min(new_list));
	return index;


def associate_domains_on_genes(bed_vals, coords, chrom, file):
	"""
	bed_vals[chrom] is a list of islands
	coords[chrom] is a list of UCSC object each of which represents a gene
	
	This is associate an island to a gene by calculating the distance between the 
	mid-point of the island and the transcription start of the gene 
	
	The algorithm runs through every gene, not particularly efficient

	Returns a list of tuples, each element of the list is a tuple made of 
	(bed_graph object, UCSC gene object)
	"""
	
	if chrom in bed_vals.keys() and chrom in coords.keys():
		association_list=[];
		genes_start_list = [];
		genes_name_list = [];
		islands = bed_vals[chrom];
		genes = coords[chrom];
		for g in genes:
			if plus.match(g.strand): genes_start_list.append(g.txStart);
			elif minus.match(g.strand): genes_start_list.append(g.txEnd);
			genes_name_list.append(g.name);
		for island_index in range(0, len(islands)):
			location = (islands[island_index].start - islands[island_index].end)/2;
			gene_index = associate_a_location_to_genes(location, genes_start_list);
			association_list.append((islands[island_index], genes[gene_index]));
			outline = islands[island_index].chrom + "\t" + genes_name_list[gene_index] + "\t" + str(islands[island_index].start) + "\t" + str(islands[island_index].end) + "\t" + str(islands[island_index].value)+'\n';
			file.write(outline);
	return association_list;


def find_location_on_genes(location, gene_regions_start_list, gene_region_end_list):
	"""
	Each region is assumed to be sequencial, namely start <= end; Even if not, 
		it will see if switching the start and end works. 
	If the location is on a gene_region, return the index
	If the location is on multiple regions, return the multiple indesis 
	If not, return a vacant list.
	
	The islands list is special in that they are non-overlapping sorted.
	The list of regions does not have to be a non-overlaping. 
	
	The algorithm runs through every gene, not particularly efficient

	"""
	assert(len(gene_regions_start_list) == len(gene_region_end_list));
	association_list = [];
	for index in range(0, len(gene_regions_start_list)):
		if location>=gene_regions_start_list[index] and location <= gene_region_end_list[index]:
			association_list.append[index];
		# incase of accidental switching of the boundaries
		elif location<=gene_regions_start_list[index] and location >= gene_region_end_list[index]:
			print "the start and end positions are not in sequencial manner";
			association_list.append[index];
	return association_list;


def find_island_on_genes(bed_graph_island, gene_regions_start_list, gene_region_end_list): 
	"""
	Each region is assumed to be sequencial, namely start <= end; Even if not, 
		it will see if switching the start and end works. 
	
	If the island is on a gene_region, return the index
	If the island is on multiple regions, return the multiple indesis 
	If not, return []
	
	The islands list is special in that they are non-overlapping sorted.
	The list of regions does not have to be a non-overlaping that 	
	
	The algorithm runs through every gene, not particularly efficient
	"""
	
	start = bed_graph_island.start;
	end = bed_graph_island.end;
	assert(len(gene_regions_start_list) == len(gene_region_end_list));
	association_list = [];
	for index in range(0, len(gene_regions_start_list)):
		if start<=gene_regions_end_list[index] and end >= gene_region_start_list[index]:
			association_list.append[index];
		# incase of accidental switching of the boundaries
		elif start<=gene_regions_start_list[index] and end >= gene_region_end_list[index]:
			print "the start and end positions are not in sequencial manner";
			association_list.append[index];
	return association_list;


def find_genic_islands(bedfile, species, known_genes, promoter_extension, region_type):
	"""
		bedfile needs to be a file with elements of BED_GRAPH type.
		known_genes is in UCSC format
		return: a BED object with elements of of BED_GRAPH type.
	"""
	if species in species_chroms.keys():
		chroms = species_chroms[species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter'];
	if region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);

	total_islands = 0;
	total_genic_islands = 0;
	bed_vals = BED.BED(species, bedfile, "BED_GRAPH", 0); #islands 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as stored one";
			sys.exit(1);
			
	coords = UCSC.KnownGenes(known_genes);
	genic_islands = {}; # a BED object of BED_GRAPH elements

	for chrom in chroms:
		genic_islands[chrom] = [];
		genes = coords[chrom];
		bed_graph_list = bed_vals[chrom];
		total_islands += len(bed_graph_list);
		# a list of flags of 0 and 1 each representing whether the corresponding island is on a gene or not
		overlapped_islands = find_domains_not_on_genes(bed_graph_list, genes, region_type, promoter_extension);
		total_genic_islands += len([elem for elem in overlapped_islands if elem != 0]);
		for index in xrange(len(overlapped_islands)):
			if overlapped_islands[index] != 0: genic_islands[chrom].append(bed_graph_list[index]);
	
	print total_genic_islands, "out of ", total_islands, " islands are inside of the " + region_type + " region";
	return genic_islands;


def place_a_gene_on_islands(island_start_list, island_end_list, region_start, region_end):
	"""
	Each region is assumed to be sequencial, namely start <= end; However, when a region is 
		not sequential, the result is 0. 
	The island list is assumed to be non-overlapping sequencial and sorted. 
	
	Returns the list of indices of whom the islands overlaps with the region.
	
	This algorithm is fast as binary search is used on island to identify overlaps.
	"""
	assert (len(island_start_list) == len(island_end_list));
	association_list = [];
	if region_start <= region_end:
		start_index = bisect.bisect_left(island_end_list, region_start);
		end_index = bisect.bisect_left(island_start_list, region_end); 
		if end_index >= start_index:
			# if start_index == len(island_end_list) the loop won't produce anything.
			for i in range(start_index, min(end_index, len(island_end_list))):	
				association_list.append(i);
	else:
		print "start < end ! not right";
	return association_list;


def gene_on_islands(island_start_list, island_end_list, region_start, region_end):
	assert (len(island_start_list) == len(island_end_list));
	assert region_start <= region_end
	start_index = bisect.bisect_left(island_end_list, region_start);
	end_index = bisect.bisect_left(island_start_list, region_end); 
	inside = 0
	if end_index > start_index:
		for i in range(start_index, min(end_index, len(island_end_list))):
			overlap_start = max(region_start, island_start_list[i])
			overlap_end = min(region_end, island_end_list[i])
			overlap_length = overlap_end - overlap_start + 1
			if overlap_length >= (region_end - region_start + 1) / 2: 
				inside = 1
	else:
		inside = -1
	return inside


def place_a_location_on_islands(island_start_list, island_end_list, position):
	"""
	Returns the list of index of whom overlaps with the position.
	This algorithm is fast as binary search is used on island to identify overlaps.
	"""
	return place_a_gene_on_islands(island_start_list, island_end_list, position, position);


def identify_gene_on_islands(bed_graph_list, genes, region_type, promoter_extension):
	#occupied_gene={};
	inside_island = []
	overlap = []
	other = []
	if len(bed_graph_list) != 0:
		island_start_list=[];
		island_end_list=[];
		#island_value_list =[];
		for item in bed_graph_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
			#island_value_list.append(item.value);
		island_start_list.sort()
		island_end_list.sort()
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
			tag = 0.0;
			indicater = gene_on_islands(island_start_list, island_end_list, start, end);
			if indicater == 1:
				inside_island.append(g.name)
			elif indicater == 0:
				overlap.append(g)
			else:
				other.append(g)
			#for i in association_list: tag += island_value_list[i];
			#occupied_gene[g.name] = tag;
	#else:
		#for g in genes: occupied_gene[g.name] = 0.0;
	return (inside_island, overlap, other);


def place_domains_on_genes(bed_graph_list, genes, region_type, promoter_extension):
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			#if g.name == '207480_s_at':
				 ##print start, end, association_list;
				 ##print island_start_list[18], island_end_list[18];
				 #print island_start_list[19], island_end_list[19];
				 #print island_start_list[20], island_end_list[20];
				 #print island_start_list[21], island_end_list[21];
			for i in association_list: tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domain_avg_score_on_genes(bed_graph_list, genes, region_type, promoter_extension):
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			if len(association_list) > 0: 
				for i in association_list:
					tag += island_value_list[i];
				occupied_gene[g.name] = tag / len(association_list);
			else:
				occupied_gene[g.name] = tag
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domains_maxscore_on_genes(bed_graph_list, genes, region_type, promoter_extension):
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			#if g.name == '207480_s_at':
				 ##print start, end, association_list;
				 ##print island_start_list[18], island_end_list[18];
				 #print island_start_list[19], island_end_list[19];
				 #print island_start_list[20], island_end_list[20];
				 #print island_start_list[21], island_end_list[21];
			for i in association_list:
				if island_value_list[i] > tag:
					tag = island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domains_maxscore_w_strand_on_genes(bed_list, genes, region_type, promoter_extension):
	"""
	region_type:
	*Promoter: TSS promoter region
	*GeneBody:  GB gene body sbtract the part already in TSS,  mutually exclusive with TSS
	*GenePromoter: G GB+TSS
	
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of gene name: tag_number
	"""
	occupied_gene={};
	if len(bed_list) != 0:
		island_start_list=[]
		island_end_list=[]
		island_value_list =[]
		island_strand_list = []
		for item in bed_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
			island_value_list.append(item.score);
			island_strand_list.append(item.strand);
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			#if g.name == '207480_s_at':
				 ##print start, end, association_list;
				 ##print island_start_list[18], island_end_list[18];
				 #print island_start_list[19], island_end_list[19];
				 #print island_start_list[20], island_end_list[20];
				 #print island_start_list[21], island_end_list[21];
			final_strand = "NA"
			for i in association_list:
				if island_value_list[i] > tag:
					tag = island_value_list[i]
					final_strand = island_strand_list[i]
			occupied_gene[g.name] = final_strand
	else:
		for g in genes: occupied_gene[g.name] = "NA";
	return occupied_gene;


def place_domain_names_on_genes(bed_graph_list, genes, region_type, promoter_extension):
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
			tag = [];
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			#if g.name == '207480_s_at':
				 ##print start, end, association_list;
				 ##print island_start_list[18], island_end_list[18];
				 #print island_start_list[19], island_end_list[19];
				 #print island_start_list[20], island_end_list[20];
				 #print island_start_list[21], island_end_list[21];
			for i in association_list:
				tag.append(island_value_list[i])
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = [];
	return occupied_gene;


def find_domainread_density_on_genes(bed_graph_list, genes, region_type, promoter_extension):
	"""
	region_type:
	*Promoter: TSS promoter region
	*GeneBody:  GB gene body sbtract the part already in TSS,  mutually exclusive with TSS
	*GenePromoter: G GB+TSS
	
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of gene_name: read_density
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list: tag += island_value_list[i];
			length = 1.0 * abs(end-start)
			occupied_gene[g.name] = tag/length;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;





def place_domain_windows_on_genes(bed_graph_list, genes, region_type, promoter_extension, threshold):
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
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list: 
				if island_value_list[i] > threshold:
					tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domain_summit_on_genes(bed_graph_list, genes, extension):
	"""
	extension: extension in bp from both ends that take count the island in
	
	returns a dictionary of gene name: summit position, which is distance from the TSS, + for downstream and - for upstream
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
			start = g.txStart - extension
			end = g.txEnd + extension
			#tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			if len(association_list) == 0:
				occupied_gene[g.name] = -1000000
			else:
				j = association_list[0]
				v = island_value_list[association_list[0]]
				for i in association_list: 
					if island_value_list[i] > v:
						j = i
						v = island_value_list[i]
				position = (island_start_list[j]+island_end_list[j]+1)/2
				#gene_length = g.txEnd - g.txStart + 1
				if plus.match(g.strand):
						occupied_gene[g.name] = position - g.txStart
				elif minus.match(g.strand):
						occupied_gene[g.name] = g.txEnd - position
	else:
		for g in genes: occupied_gene[g.name] = -1000000;
	return occupied_gene;


def place_domains_on_promoters(bed_graph_list, genes, start_extension, end_extension):
	"""
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of region name: tag_number
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
			if plus.match(g.strand):
				start = g.txStart - start_extension
				end = g.txStart + end_extension
			# When strand is negative, the regions are still sequential
			elif minus.match(g.strand):
				start = g.txEnd - end_extension
				end = g.txEnd + start_extension
			else: print "Gene strand identification crisis!";
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list:
				tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domains_on_gene_ends(bed_graph_list, genes, start_extension, end_extension):
	"""
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of region name: tag_number
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
			if plus.match(g.strand):
				start = g.txEnd - start_extension
				end = g.txEnd + end_extension
			# When strand is negative, the regions are still sequential
			elif minus.match(g.strand):
				start = g.txStart - end_extension
				end = g.txStart + start_extension
			else: print "Gene strand identification crisis!";
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list:
				tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


### didn't consider strands!!!
def place_domains_on_regions(bed_graph_list, genes, start_extension, end_extension, threshold):
	"""
	returns a dictionary of region name: tag_number
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
			start = g.txStart - start_extension
			end = g.txEnd + end_extension
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list:
				if island_value_list[i] > threshold:
					tag += island_value_list[i];
			occupied_gene[g.name] = tag;
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_domains_on_regions_with_changing_thresholds(bed_graph_list, genes, start_extension, end_extension, windowsize, thresholdtable):
	"""
	promoter_extension: the value for the upstream extension
	
	returns a dictionary of region name: tag_number
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
			start = g.txStart - start_extension
			end = g.txEnd + end_extension
			window = (end - start)/windowsize + 2
			threshold = thresholdtable[window]
			tag = 0.0;
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list:
				tag += island_value_list[i]
			if tag > threshold:
				occupied_gene[g.name] = tag;
			else:
				occupied_gene[g.name] = 0.0
	else:
		for g in genes: occupied_gene[g.name] = 0.0;
	return occupied_gene;


def place_summary_on_enhancers(bed_graph_list, genes, start_extension, end_extension, windowsize):
	result_value_list = []
	if len(bed_graph_list) != 0:
		island_start_list=[];
		island_end_list=[];
		island_value_list =[];
		for item in bed_graph_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
			island_value_list.append(item.value);
		for g in genes:
			start = g.txStart - start_extension
			end = g.txEnd + end_extension
			association_list = place_a_gene_on_islands(island_start_list, island_end_list, start, end);
			for i in association_list:
				result_value_list.append(island_value_list[i])
	return result_value_list


def output_marked_geneset_in_UCSC(species, bedfile, known_genes_file, region_type, promoter_extension, outfile):
	"""
	Output a file of marked genes in UCSC format.
	Input: 
		bedfile: bed file with bed_graph elements, such as islands. 
		know_gene_file: UCSC type file 
		region type: 'Promoter', 'GeneBody', 'GenePromoter'
	"""
	if species in species_chroms.keys():
		chroms = species_chroms[species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter'];
	if region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);
		
	bed_vals = BED.BED(species, bedfile, "BED_GRAPH", 0); #islands 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(known_genes_file);
	occupied_genes = {}; 
	total = 0.0;

	for chrom in chroms:
		genes = coords[chrom];
		total += len(genes);
		bed_graph_list = bed_vals[chrom];
		current_genes = place_domains_on_genes(bed_graph_list, genes, region_type, promoter_extension);
		for gene in current_genes.keys():
			if current_genes[gene] >.01:
				occupied_genes[gene] = current_genes[gene];
	print "Number of marked genes is: ", len(occupied_genes), " out of total of ", total, " genes in "+ species ;
	
	gene_list = occupied_genes.keys();
	gene_set_manipulation.output_UCSCsubset_in_file(known_genes_file, gene_list, outfile);


def find_islands_not_on_genenames(island_start_list, island_end_list, gene_start_list, gene_end_list, gene_name_list):
	"""
	Returns the list of island each of which has a flag telling whether it overlaps our not.
	This algorithm is fast as binary search is used on island to identify overlaps.
	"""
	assert (len(island_start_list) == len(island_end_list));
	assert (len(gene_start_list) == len(gene_end_list));
	overlapped_islands=[0]*len(island_start_list);

	for gene_index in range(0, len(gene_start_list)):
		region_start = gene_start_list[gene_index];
		region_end = gene_end_list[gene_index];
		association_list = place_a_gene_on_islands(island_start_list, island_end_list, region_start, region_end);
		#print gene_name_list[gene_index], association_list;
		#for item in association_list:	
		#	print gene_name_list[gene_index], island_start_list[item], island_end_list[item];
		for item in association_list:
			overlapped_islands[item] = 1;
	return 	overlapped_islands


def find_islands_not_on_genes(island_start_list, island_end_list, gene_start_list, gene_end_list):
	"""
	Returns the list of island each of which has a flag telling whether it overlaps our not.
	This algorithm is fast as binary search is used on island to identify overlaps.
	"""
	assert (len(island_start_list) == len(island_end_list));
	assert (len(gene_start_list) == len(gene_end_list));
	overlapped_islands=[0]*len(island_start_list);

	for gene_index in range(0, len(gene_start_list)):
		region_start = gene_start_list[gene_index];
		region_end = gene_end_list[gene_index];
		association_list = place_a_gene_on_islands(island_start_list, island_end_list, region_start, region_end);
		#print gene_index, association_list;
		for item in association_list:
			overlapped_islands[item] = 1;
	return 	overlapped_islands


def find_domains_not_on_genes(bed_graph_list, genes, region_type, promoter_extension):
	"""
	region_type:
	*Promoter: TSS promoter region
	*GeneBody:  GB gene body sbtract the part already in TSS,  mutually exclusive with TSS
	*GenePromoter: G GB+TSS
	
	promoter_extension: the value for the upstream extension
	
	returns the list of islands each with a flag, if o it is not on a gene.
	"""
	overlapped_islands={};
	if len(bed_graph_list) != 0:
		island_start_list=[];
		island_end_list=[];
		for item in bed_graph_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
		region_start_list = [];
		region_end_list = [];
		gene_name_list=[];
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
			region_start_list.append(start);
			region_end_list.append(end);
			gene_name_list.append(g.name);
		overlapped_islands = find_islands_not_on_genenames(island_start_list, island_end_list, region_start_list, region_end_list, gene_name_list);	
	return overlapped_islands;	


def output_bed_graph(bed_vals, file_name, species):
	"""
	print out a bed_val dictionary made of bed_graph objects.
	"""
	file = open(file_name, 'w');
	chroms = species_chroms[species];
	chroms.sort();
	for chrom in chroms:
		for bed_graph in bed_vals[chrom]:
			outline = bed_graph.chrom + ' ' + str(bed_graph.start) + ' ' + str(bed_graph.end) + ' ' + str(bed_graph.value)+'\n';
			file.write(outline);
	file.close();


def find_integenic_islands(bedfile, species, known_genes, promoter_extension, region_type):
	"""
		bedfile needs to be a file with elements of BED_GRAPH type.
		known_genes is in UCSC format
		return: a BED object with elements of of BED_GRAPH type.
	"""
	if species in species_chroms.keys():
		chroms = species_chroms[species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter'];
	if region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);

	total_islands = 0;
	total_integenic_islands = 0;
	bed_vals = BED.BED(species, bedfile, "BED_GRAPH", 0); #islands 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as stored one";
			sys.exit(1);
			
	coords = UCSC.KnownGenes(known_genes);
	integenic_islands = {}; # a BED object of BED_GRAPH elements

	for chrom in chroms:
		integenic_islands[chrom] = [];
		genes = coords[chrom];
		bed_graph_list = bed_vals[chrom];
		total_islands += len(bed_graph_list);
		# a list of flags of 0 and 1 each representing whether the corresponding island is on a gene or not
		overlapped_islands = find_domains_not_on_genes(bed_graph_list, genes, region_type, promoter_extension);
		total_integenic_islands += len([elem for elem in overlapped_islands if elem==0]);
		for index in xrange(len(overlapped_islands)):
			if overlapped_islands[index] == 0: integenic_islands[chrom].append(bed_graph_list[index]);
	
	print total_integenic_islands, "out of ", total_islands, " islands are outside of the " + region_type + " region";
	return integenic_islands;


def test1():
	# module testing
	island_start_list= [1,3, 15, 21, 37, 111, 200];
	island_end_list = [2, 5, 20, 30, 100, 151, 240];
	gene_start_list = [1, 30, 400];
	gene_end_list = [10, 150, 450];
	overlapped_islands = find_islands_not_on_genes(island_start_list, island_end_list, gene_start_list, gene_end_list);
	for i in xrange(len(overlapped_islands)):
		if overlapped_islands[i] == 0: print island_start_list[i], island_end_list[i];
		
	print associate_a_location_to_genes(111, island_start_list)

def get_gene_bounded_region_lists(genes, chrom, bound_file, extension):
	TSS_list = []
	region_start_list = []
	region_end_list = []
	gene_name_list = []
	start_list = []
	end_list = []
	f = open(bound_file,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			start_list.append(atoi(sline[1]))
			end_list.append(atoi(sline[2]))
	f.close()
	start_list.sort()
	end_list.sort()
	for g in genes:
		if plus.match(g.strand):
			TSS = g.txStart
		elif minus.match(g.strand):
			TSS = g.txEnd
		else: print "Gene strand identification crisis!"
		start_index = bisect.bisect_left(start_list, TSS)
		end_index = bisect.bisect_left(end_list, TSS)
		if start_index - end_index == 1:
			start = max(TSS - extension, 0)
			end = TSS + extension
		elif start_index == end_index and start_index > 0 and start_index < len(start_list):
			start = max(end_list[end_index - 1], TSS - extension)
			end = min(start_list[start_index], TSS + extension)
		elif start_index == 0 and len(start_list) > 0:
			start = max(TSS - extension, 0)
			end = min(start_list[start_index], TSS + extension)
		elif start_index == len(start_list) and len(start_list) > 0:
			start = max(end_list[end_index - 1], TSS - extension)
			end = TSS + extension
		else: 
			start = max(TSS - extension, 0)
			end = TSS + extension
		TSS_list.append(TSS)
		region_start_list.append(start)
		region_end_list.append(end)
		gene_name_list.append(g.name)
	return gene_name_list, TSS_list, region_start_list, region_end_list


def get_gene_feature_region_lists(genes, region_type, promoter_extension):
	region_start_list = [];
	region_end_list = [];
	gene_name_list=[];
	for g in genes:
		# gene body and TSS are mutually exclusive
		if plus.match(g.strand):
			TSS_start = g.txStart - promoter_extension
			TSS_end = g.txStart + promoter_extension
			GB_start = TSS_end
			GB_end = g.txEnd
			G_start = TSS_start
			G_end = GB_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension
		# When strand is negative, the regions are still sequential
		elif minus.match(g.strand):
			TSS_start = g.txEnd - promoter_extension
			TSS_end = g.txEnd + promoter_extension
			GB_start = g.txStart
			GB_end = TSS_start
			G_start = GB_start
			G_end = TSS_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension
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
		elif region_type == 'GeneExtended':
			start = GE_start
			end = GE_end
		region_start_list.append(start);
		region_end_list.append(end);
		gene_name_list.append(g.name);
	return gene_name_list, region_start_list, region_end_list


def get_promoter_feature_region_lists(genes, start_extension, end_extension):
	region_start_list = [];
	region_end_list = [];
	gene_name_list=[];
	for g in genes:
		# gene body and TSS are mutually exclusive
		if plus.match(g.strand):
			TSS_start = g.txStart - start_extension
			TSS_end = g.txStart + end_extension
			'''GB_start = TSS_end
			GB_end = g.txEnd
			G_start = TSS_start
			G_end = GB_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension'''
		# When strand is negative, the regions are still sequential
		elif minus.match(g.strand):
			TSS_start = g.txEnd - end_extension
			TSS_end = g.txEnd + start_extension
			'''GB_start = g.txStart
			GB_end = TSS_start
			G_start = GB_start
			G_end = TSS_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension'''
		else: print "Gene strand identification crisis!";
		
		'''if region_type == 'Promoter':
			start = TSS_start
			end = TSS_end
		elif region_type == 'GeneBody': # start could be bigger than end if TSS is chosen big.
			start = GB_start
			end = GB_end
		elif region_type == 'GenePromoter':
			start = G_start
			end = G_end
		elif region_type == 'GeneExtended':
			start = GE_start
			end = GE_end'''
		region_start_list.append(TSS_start);
		region_end_list.append(TSS_end);
		gene_name_list.append(g.name);
	return gene_name_list, region_start_list, region_end_list


def get_geneend_feature_region_lists(genes, start_extension, end_extension):
	region_start_list = [];
	region_end_list = [];
	gene_name_list=[];
	for g in genes:
		# gene body and TSS are mutually exclusive
		if plus.match(g.strand):
			TSS_start = g.txEnd - start_extension
			TSS_end = g.txEnd + end_extension
			'''GB_start = TSS_end
			GB_end = g.txEnd
			G_start = TSS_start
			G_end = GB_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension'''
		# When strand is negative, the regions are still sequential
		elif minus.match(g.strand):
			TSS_start = g.txStart - end_extension
			TSS_end = g.txStart + start_extension
			'''GB_start = g.txStart
			GB_end = TSS_start
			G_start = GB_start
			G_end = TSS_end
			GE_start = g.txStart - promoter_extension
			GE_end = g.txEnd + promoter_extension'''
		else: print "Gene strand identification crisis!";
		
		'''if region_type == 'Promoter':
			start = TSS_start
			end = TSS_end
		elif region_type == 'GeneBody': # start could be bigger than end if TSS is chosen big.
			start = GB_start
			end = GB_end
		elif region_type == 'GenePromoter':
			start = G_start
			end = G_end
		elif region_type == 'GeneExtended':
			start = GE_start
			end = GE_end'''
		region_start_list.append(TSS_start);
		region_end_list.append(TSS_end);
		gene_name_list.append(g.name);
	return gene_name_list, region_start_list, region_end_list


def find_closest_island_around_promoter(bed_graph_list, genes):
	"""
	promoter_extension: the value for the upstream extension
	returns a dictionary of gene name: tag_number
	islands must be sorted!
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
		#island_start_list.sort()
		#island_end_list.sort()
		for g in genes:
			if plus.match(g.strand):
				TSS = g.txStart
				start_index = bisect.bisect_left(island_start_list, TSS)
				end_index = bisect.bisect_left(island_end_list, TSS)
				if start_index - end_index == 1:  # TSS in an island
					occupied_gene[g.name] = "0\t" + str(island_value_list[end_index]) + "\t0\t" + str(island_value_list[end_index]) + "\t" + str(abs(g.txStart - g.txEnd));
				else:
					assert end_index == start_index
					#if end_index != start_index:
					#	print g.name, str(end_index - start_index);
					index = end_index
					if index == 0:
						down_distance = island_start_list[index] - TSS
						assert down_distance >= 0
						down_value = island_value_list[index]
						occupied_gene[g.name] = "None\tNone\t" + str(down_distance) + "\t" + str(down_value) + "\t" + str(abs(g.txStart - g.txEnd));
					elif index == len(island_end_list):
						index = index - 1
						up_distance = TSS - island_end_list[index]
						assert up_distance >= 0
						up_value = island_value_list[index]
						occupied_gene[g.name] =  str(up_distance) + "\t" + str(up_value) + "\tNone\tNone\t" + str(abs(g.txStart - g.txEnd));
					else:
						up_index = index - 1
						up_distance = TSS - island_end_list[up_index]
						up_value = island_value_list[up_index]
						down_distance = island_start_list[index] - TSS
						down_value = island_value_list[index]
						assert up_distance >= 0
						assert down_distance >= 0
						occupied_gene[g.name] =  str(up_distance) + "\t" + str(up_value) + "\t" + str(down_distance) + "\t" + str(down_value) + "\t" + str(abs(g.txStart - g.txEnd));
			elif minus.match(g.strand):
				TSS = g.txEnd
				start_index = bisect.bisect_left(island_start_list, TSS)
				end_index = bisect.bisect_left(island_end_list, TSS)
				if start_index - end_index == 1:  # TSS in an island
					occupied_gene[g.name] = "0\t" + str(island_value_list[end_index]) + "\t0\t" + str(island_value_list[end_index]) + "\t" + str(abs(g.txStart - g.txEnd));
				else:
					assert end_index == start_index
					#if end_index != start_index:
					#	print g.name, str(end_index - start_index);
					index = end_index
					if index == 0:
						up_distance = island_start_list[index] - TSS
						assert up_distance >= 0
						up_value = island_value_list[index]
						occupied_gene[g.name] = str(up_distance) + "\t" + str(up_value) + "\tNone\tNone\t" + str(abs(g.txStart - g.txEnd));
					elif index == len(island_end_list):
						index = index - 1
						down_distance = TSS - island_end_list[index]
						assert down_distance >= 0
						down_value = island_value_list[index]
						occupied_gene[g.name] =  "None\tNone\t" + str(down_distance) + "\t" + str(down_value) + "\t" + str(abs(g.txStart - g.txEnd));
					else:
						down_index = index - 1
						down_distance = TSS - island_end_list[down_index]
						down_value = island_value_list[down_index]
						up_distance = island_start_list[index] - TSS
						up_value = island_value_list[index]
						assert up_distance >= 0
						assert down_distance >= 0
						occupied_gene[g.name] =  str(up_distance) + "\t" + str(up_value) + "\t" + str(down_distance) + "\t" + str(down_value) + "\t" + str(abs(g.txStart - g.txEnd));
			else: print "Gene strand identification crisis!";
	else:
		for g in genes: occupied_gene[g.name] = "None\tNone\tNone\tNone\t" + str(abs(g.txStart - g.txEnd));
	return occupied_gene;


def find_closest_island_distance_around_gene(bed_graph_list, genes):
	"""
	promoter_extension: the value for the upstream extension
	returns a dictionary of gene name: tag_number
	islands must be sorted!
	"""
	occupied_gene={};
	if len(bed_graph_list) != 0:
		island_start_list=[];
		island_end_list=[];
		#island_value_list =[];
		for item in bed_graph_list:
			island_start_list.append(item.start);
			island_end_list.append(item.end);
		island_start_list.sort()
		island_end_list.sort()
			#island_value_list.append(item.value);
		for g in genes:
			if plus.match(g.strand):
				TSS = g.txStart
				start_index = bisect.bisect_left(island_start_list, TSS)
				end_index = bisect.bisect_left(island_end_list, TSS)
				if start_index - end_index == 1:  # TSS in an island
					occupied_gene[g.name] = "0";
				else:
					assert end_index == start_index
					#if end_index != start_index:
					#	print g.name, str(end_index - start_index);
					index = end_index
					if index == 0:
						down_distance = island_start_list[index] - TSS
						assert down_distance >= 0
						gene_length = abs(g.txStart - g.txEnd)
						distance = down_distance - gene_length
						if distance < 0:
							distance = 0
						#own_value = island_value_list[index]
						occupied_gene[g.name] = str(distance);
					elif index == len(island_end_list):
						index = index - 1
						up_distance = TSS - island_end_list[index]
						assert up_distance >= 0
						#up_value = island_value_list[index]
						occupied_gene[g.name] =  str(up_distance);
					else:
						up_index = index - 1
						up_distance = TSS - island_end_list[up_index]
						#up_value = island_value_list[up_index]
						down_distance = island_start_list[index] - TSS
						#down_value = island_value_list[index]
						gene_length = abs(g.txStart - g.txEnd)
						down_distance = down_distance - gene_length
						if down_distance < 0:
							down_distance = 0
						assert up_distance >= 0
						assert down_distance >= 0
						occupied_gene[g.name] =  str(min(up_distance, down_distance));
			elif minus.match(g.strand):
				TSS = g.txEnd
				start_index = bisect.bisect_left(island_start_list, TSS)
				end_index = bisect.bisect_left(island_end_list, TSS)
				if start_index - end_index == 1:  # TSS in an island
					occupied_gene[g.name] = "0\t" + str(island_value_list[end_index]) + "\t0\t" + str(island_value_list[end_index]) + "\t" + str(abs(g.txStart - g.txEnd));
				else:
					assert end_index == start_index
					#if end_index != start_index:
					#	print g.name, str(end_index - start_index);
					index = end_index
					if index == 0:
						up_distance = island_start_list[index] - TSS
						assert up_distance >= 0
						#up_value = island_value_list[index]
						occupied_gene[g.name] = str(up_distance);
					elif index == len(island_end_list):
						index = index - 1
						down_distance = g.txStart - island_end_list[index]
						if down_distance < 0:
							down_distance = 0
						assert down_distance >= 0
						#down_value = island_value_list[index]
						occupied_gene[g.name] =  str(down_distance);
					else:
						down_index = index - 1
						down_distance = g.txStart - island_end_list[down_index]
						if down_distance < 0:
							down_distance = 0
						#down_value = island_value_list[down_index]
						up_distance = island_start_list[index] - TSS
						#up_value = island_value_list[index]
						assert up_distance >= 0
						assert down_distance >= 0
						occupied_gene[g.name] =  str(min(up_distance, down_distance));
			else: print "Gene strand identification crisis!";
	else:
		for g in genes: occupied_gene[g.name] = "None";
	return occupied_gene;


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="islands file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in")
	parser.add_option("-f", "--outfileoccupiedgenes", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for occupied genes")
	
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
		
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0); #islands 
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; # a BED object of BED_GRAPH elements
	total = 0;

	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			total += len(genes);
			bed_graph_list = bed_vals[chrom];
			current_genes = place_domains_on_genes(bed_graph_list, genes, opt.region_type, opt.promoter_extension);
			for gene in current_genes.keys():
				if current_genes[gene] >.01:
					occupied_genes[gene] = current_genes[gene];
	print "Number of marked genes is: ", len(occupied_genes), " out of total of ", total, " genes in "+ opt.species ;
	
	gene_list = occupied_genes.keys();
	gene_set_manipulation.output_UCSCsubset_in_file(opt.known_genes, gene_list, opt.out_file);


if __name__ == "__main__":
	main(sys.argv)

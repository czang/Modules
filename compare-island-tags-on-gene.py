#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;

import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is used to associate the islands distribution with genes comparing between two cell types. For each known gene along a genome, it tells the ratio of tag counts on islands between the two cell types. 

The creteria that an island is in a gene region is that the mean position (defined as the center of mass) of the island is in the gene region.
"""


def getTagCounts(bed_vals, chrom):
        """
		bed should be a bed object with the bed_graph data type
		The data structure is a dictionary keyed by chromosome and valued by   
				a bunch of bed_graph object each of which contains [chrom start end count]
        Return a list of TagCounts on a given chromosome
        """
        TagCounts = [];
        for t in bed_vals[chrom]:
            TagCounts.append(t.value);
        try:
            return TagCounts;
        except:
            sys.stderr.write("Having trouble returning ends %s\n" % TagCounts)
            return ''


def indexsearch(list, search):    
    """ classic binary search -- If 'search' is in 'list', its index is
    returned.  If not, the index where 'search' should be in 'list' is returned """
    right = len(list)
    left = 0
    previous_center = 0
    if search < list[0]:
        return 0
    while 1:
        center = (left + right) / 2
        candidate = list[center]
        if search == candidate:
            return center
        if center == previous_center:
            return (1 + center);
        elif search < candidate:
            right = center
        else:
            left = center
        previous_center = center


def FindIslandCenter(IslandPosition, StartList, EndList, TagCountList, WindowSize):
	'''Return the average central coordinate of an island, with window counts weighted'''
	StartPosition = indexsearch(StartList, IslandPosition['start'])
	EndPosition = indexsearch(EndList, IslandPosition['end'])
	average = 0.0
	Total = 0.0;
	TotalCounts = 0.0;
	for i in range(StartPosition, EndPosition+1):
		Total = Total + TagCountList[i] * ( StartList[i] + WindowSize / 2 )
		TotalCounts = TotalCounts + TagCountList[i]
	average = Total / TotalCounts
	return int(average)


def get_total_tag_count_on_islands(islands):
	total_counts = 0
	for t in islands:
		total_counts += t.value
	return total_counts

def chrom_tag_density(Histogram_bed, chrom, length):
	tag_list = getTagCounts(Histogram_bed, chrom)
	total_tag = 0
	for tag in tag_list:
		total_tag += tag
	return float(total_tag) / float(length)
	

def get_island_mean_position_list(bed_vals, Histogram_bed, chrom, window_size):
	'''Return a list of mean positions of islands in one chromosome'''
	island_mean_position_list = []
	if chrom in bed_vals.keys():
		StartList = Histogram_bed.getStarts(chrom)
		EndList = Histogram_bed.getEnds(chrom)
		TagCountList = getTagCounts(Histogram_bed, chrom)
		islands = bed_vals[chrom]
		IslandPosition = {}
		for island in islands: 
			IslandPosition['start'] = island.start
			IslandPosition['end'] = island.end
			island_mean_position_list.append(FindIslandCenter(IslandPosition, StartList, EndList, TagCountList, window_size))
	return island_mean_position_list


def get_tag_number_in_genes(island_position_list, island_value_list, genes, TSSthreshold, region, psuedocount):
	gene_tag_number = {}
	if len(island_position_list) != 0:
		for g in genes:
			if g.strand == '+':
				TSS_start = g.txStart - TSSthreshold
				TSS_end = g.txStart + TSSthreshold
				GB_start = TSS_end
				GB_end = g.txEnd
				G_start = TSS_start
				G_end = GB_end
			elif g.strand == '-':
				TSS_start = g.txEnd - TSSthreshold
				TSS_end = g.txEnd + TSSthreshold
				GB_start = g.txStart
				GB_end = TSS_start
				G_start = GB_start
				G_end = TSS_end
			if region == 'TSS':
				start = TSS_start
				end = TSS_end
			elif region == 'GeneBody':
				start = GB_start
				end = GB_end
			elif region == 'Gene':
				start = G_start
				end = G_end
			start_index = indexsearch(island_position_list, start)
			end_index = indexsearch(island_position_list, end)
			assert (end_index >= start_index)
			if end_index == start_index:
				gene_tag_number[g.name] = psuedocount
			else: 
				tag = 0
				for i in range(start_index, end_index):
					tag = tag + island_value_list[i]
				gene_tag_number[g.name] = tag
	else:
		for g in genes:
			gene_tag_number[g.name] = psuedocount
	return gene_tag_number


def result_comparison(Result_1, Result_2, psuedocount, normalization):
	result = {}
	if len(Result_1.keys()) != len(Result_2.keys()):
		print 'result error'
	for name in Result_1.keys():
		if Result_2.has_key(name):
			if Result_1[name] == psuedocount:
				if Result_2[name] != psuedocount:
					Result_2[name] = Result_2[name] + psuedocount
			else:
				Result_1[name] = float(Result_1[name]) / float(normalization)
				if Result_2[name] == psuedocount:
					Result_1[name] = Result_1[name] + psuedocount
			result[name] = log(float(Result_1[name]) / float(Result_2[name]),2)
	return result



def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="file1 with island info in BED format")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="file2 with island info in BED format")
	parser.add_option("-d", "--bedfile3", action="store", type="string", dest="bedfile3", metavar="<file>", help="file1 with summary info in bed graph format")
	parser.add_option("-g", "--bedfile4", action="store", type="string", dest="bedfile4", metavar="<file>", help="file2 with summary info in bed graph format")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-w", "--window_size", action="store", type="int", dest="windowsize", help="window size of summary", metavar="<int>")
	parser.add_option("-t", "--TSS_threshold", action="store", type="int", dest="TSSthreshold", help="upstream and downstream threshold of TSS region", metavar="<int>")
	parser.add_option("-p", "--pseudocount", action="store", type="int", dest="pseudocount", help="pseudo tag count for a gene that has no island in it", metavar="<int>")
	parser.add_option("-r", "--'TSS' or 'GeneBody' or 'Gene'", action="store", type="string", dest="region", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 22:
        	parser.print_help()
        	sys.exit(1)
	
	mm_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM']
	hg_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
	
	if opt.species == "mm8": 
		chroms = mm_chroms
	elif opt.species == "hg18": 
		chroms = hg_chroms
	else:
		print 'Species Error'


	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0)
	bed_vals_2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0)
	coords = UCSC.KnownGenes(opt.known_genes)
	
	total_tag_1 = 0
	total_tag_2 = 0
	for chrom in chroms: 
		total_tag_1 += get_total_tag_count_on_islands(bed_vals_1[chrom])
		total_tag_2 += get_total_tag_count_on_islands(bed_vals_2[chrom])
	Normalization = float(total_tag_1) / float(total_tag_2)
	
	Result = {}
	
	for chrom in chroms: 
		Histogram_bed_1 = BED.BED(opt.species, opt.bedfile3 + '_temp_' + chrom, "BED_GRAPH")
		Histogram_bed_2 = BED.BED(opt.species, opt.bedfile4 + '_temp_' + chrom, "BED_GRAPH")
		island_list_1 = get_island_mean_position_list(bed_vals_1, Histogram_bed_1, chrom, opt.windowsize)
		island_list_2 = get_island_mean_position_list(bed_vals_2, Histogram_bed_2, chrom, opt.windowsize)
		Result_1 = {}
		Result_2 = {}
		if chrom in coords.keys():
			Result_1 = get_tag_number_in_genes(island_list_1, getTagCounts(bed_vals_1, chrom), coords[chrom], opt.TSSthreshold, opt.region, opt.pseudocount)
		if chrom in coords.keys():
			Result_2 = get_tag_number_in_genes(island_list_2, getTagCounts(bed_vals_2, chrom), coords[chrom], opt.TSSthreshold, opt.region, opt.pseudocount)
		Result[chrom] = result_comparison(Result_1, Result_2, opt.pseudocount, Normalization)
	
	f = open (opt.outfile,'w')
	for chrom in Result.keys():
		for g in Result[chrom].keys():
			f.write(g + '\t' + str(Result[chrom][g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
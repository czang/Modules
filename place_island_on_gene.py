#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;

import UCSC;
from species_chroms import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is used to associate the islands distribution with genes. 
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
	'''
	Return the average central coordinate of an island, with window counts weighted

	StartPosition: the first tagged window in the island
	EndPosition: the last tagged window in the island

	For a regular island, the start(end) position shall be the start of a window in the StartList(EndList)
	In general, it might not be true. The bed summary file saves the sorted windows with tags in them,
	Two windows might be far away from each other.
	'''
	StartPosition = indexsearch(StartList, IslandPosition['start'])
	if (StartPosition > 0) and (IslandPosition['start'] <=  EndList[StartPosition-1] ):
		StartPosition = StartPosition - 1;
		print "In FindIslandCenter a special situation happens";
	EndPosition = indexsearch(EndList, IslandPosition['end'])
	if (EndPosition == len(EndList)):
		EndPosition = EndPosition - 1;
		print "In FindIslandCenter a special situation happens";
	elif (IslandPosition['end'] < StartList[EndPosition]):
		EndPosition = EndPosition - 1;
		print "In FindIslandCenter a special situation happens";
	average = 0.0
	Total = 0.0;
	TotalCounts = 0.0;
	for i in range(StartPosition, EndPosition+1):
		Total = Total + TagCountList[i] * ( StartList[i] + WindowSize / 2 )
		TotalCounts = TotalCounts + TagCountList[i]
	average = float(Total) / float(TotalCounts)
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
	else: print chrom + " is not in recognized !!!";
	return island_mean_position_list



def get_tag_number_in_genes(island_position_list, island_value_list, genes, TSSthreshold, region):
	"""
	Find out for each gene, how many islands in island_position_list are in the defined region.
	"""
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
			if region == 'Promoter':
				start = TSS_start
				end = TSS_end
			elif region == 'GeneBody': # start could be bigger than end if TSS is chosen big.
				start = GB_start
				end = GB_end
			elif region == 'GenePromoter':
				start = G_start
				end = G_end
			start_index = indexsearch(island_position_list, start)
			end_index = indexsearch(island_position_list, end)
			if end_index <= start_index: 
				gene_tag_number[g.name] = 0
			else: 
				tag = 0
				for i in range(start_index, end_index):
					tag = tag + island_value_list[i]
				gene_tag_number[g.name] = tag
	else:
		for g in genes:
			gene_tag_number[g.name] = 0
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
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="file with island info in BED format")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="file with summary info in bed graph format")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-w", "--window_size", action="store", type="int", dest="windowsize", help="window size of summary", metavar="<int>")
	parser.add_option("-t", "--TSS_threshold", action="store", type="int", dest="TSSthreshold", help="upstream and downstream threshold of TSS region", metavar="<int>")
	parser.add_option("-r", "--'Promoter' or 'GeneBody' or 'GenePromoter'", action="store", type="string", dest="region", metavar="<str>", help="region to count tags in")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
        	parser.print_help()
        	sys.exit(1)
	
	if species_chroms.has_key(opt.species): 
		chroms = species_chroms[opt.species];
	else:
		print 'Species not recognized'


	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0) # islands
	coords = UCSC.KnownGenes(opt.known_genes)
	
	Result = {}
	
	for chrom in chroms: 
		Histogram_bed_1 = BED.BED(opt.species, opt.bedfile2 + '_temp_' + chrom, "BED_GRAPH")
		island_list_1 = get_island_mean_position_list(bed_vals_1, Histogram_bed_1, chrom, opt.windowsize)
		Result_1 = {}
		if chrom in coords.keys():
			Result_1 = get_tag_number_in_genes(island_list_1, getTagCounts(bed_vals_1, chrom), coords[chrom], opt.TSSthreshold, opt.region)
		Result[chrom] = Result_1
	
	f = open (opt.outfile,'w')
	for chrom in Result.keys():
		for g in Result[chrom].keys():
			if Result[chrom][g] > 1:
				f.write(g + '\t' + str(Result[chrom][g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)

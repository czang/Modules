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

## get BED module
import BED;

import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();



def getTagCounts(bed, chrom):
        """
		bed should be a bed object with the bed_graph data type
		The data structure is a dictionary keyed by chromosome and valued by   
				a bunch of bed_graph object each of which contains [chrom start end count]
        Return a list of TagCounts on a given chromosome
        """
        TagCounts = [];
        for t in bed.bed_vals[chrom]:
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


"""
This module is used to associate the islands distribution with genes. Of all islands through the genome of a cell, it tells: 
1) How many are in the promotor region 
2) How many are in the gene body region 
3) In gene body region, how many are in exons and how many are in introns, respect with the bp length of exons and introns 
4) How many are in the intergenic region

"""


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


def aver_tag_density_on_islands(islands):
	"""return the average tag density of the islands in one chromosome. It is the ratio between the total tag number in all islands and the total island lengths """
	total_bp = 0.0
	total_counts = 0.0
	for island_index in range(0, len(islands)):
		total_bp = total_bp - islands[island_index].start + islands[island_index].end + 1
		total_counts = total_counts + islands[island_index].value
	return float(total_counts) / float(total_bp)


def get_island_regions(bed_vals, Histogram_bed, coords, chrom, TSSthreshold, window_size):
	Result = {}
	Result['total'] = 0
	Result['TSS'] = 0
	Result['GeneBody'] = 0
	Result['exon'] = 0
	Result['Intergenic'] = 0
	Result['exon_length'] = 0.0
	Result['gene_length'] = 0.0

	if chrom in bed_vals.keys() and chrom in coords.keys():
		if window_size > 0:
			StartList = Histogram_bed.getStarts(chrom)
			EndList = Histogram_bed.getEnds(chrom)
			TagCountList = getTagCounts(Histogram_bed, chrom)
		islands = bed_vals[chrom]
		
		genes = coords[chrom]
		TSS_List = [] ##make TSS list
		Start_List = []
		End_List = []
		
		for g in genes: 
			Start_List.append(g.txStart)
			End_List.append(g.txEnd)
			Result['gene_length'] = Result['gene_length'] + g.txEnd - g.txStart + 1
			'''total_exon_Start_List = [] 
			total_exon_End_List = [] '''
			total_exon_Starts = g.exonStarts.replace(',',' ')
			total_exon_Ends = g.exonEnds.replace(',',' ')
			total_exon_Starts = total_exon_Starts.split()
			total_exon_Ends = total_exon_Ends.split()
			'''for item in total_exon_Starts:
				total_exon_Start_List.append(int(item))
			for item in total_exon_Ends:
				total_exon_End_List.append(int(item))'''
			for i in range(0, len(total_exon_Starts)):
				Result['exon_length'] = Result['exon_length'] + float(total_exon_Ends[i]) - float(total_exon_Starts[i]) + 1
			
			if plus.match(g.strand):
				TSS_List.append(g.txStart);
			elif minus.match(g.strand):
				TSS_List.append(g.txEnd)
		TSS_List.sort()
		
		IslandPosition = {}
		average = 0
		TSSindex = 0
		
		TSS_Number = 0
		GeneBody_Number = 0
		exon_Number = 0
		Intergenic_Number = 0
		
		for island_index in range(0, len(islands)):
			IslandPosition['start'] = islands[island_index].start
			IslandPosition['end'] = islands[island_index].end
			if window_size > 0:
				average = FindIslandCenter(IslandPosition, StartList, EndList, TagCountList, window_size)
			else:
				average = int((IslandPosition['start']+IslandPosition['end'])/2)
			##print chrom, islands[island_index].start, islands[island_index].end, average
			TSSindex = indexsearch(TSS_List, average)
			index_start = indexsearch(Start_List, average)
			index_end = indexsearch(End_List, average)
			if abs(average - TSS_List[min(TSSindex,(len(TSS_List)-1))]) < TSSthreshold or abs(average - TSS_List[max(0,(TSSindex-1))]) < TSSthreshold:
				TSS_Number += 1; 
				#if abs(average - TSS_List[min(TSSindex,(len(TSS_List)-1))]) < TSSthreshold:
					#print genes[min(TSSindex,(len(TSS_List)-1))].name
				#else:
					#print genes[max(0,(TSSindex-1))].name
			elif index_start != index_end: 
				GeneBody_Number += 1;
				#print genes[index_end].name
				exon_Start_List = []
				exon_End_List = []
				exon_Starts = genes[index_end].exonStarts.replace(',',' ')
				exon_Ends = genes[index_end].exonEnds.replace(',',' ')
				exon_Starts = exon_Starts.split()
				exon_Ends = exon_Ends.split()
				for item in exon_Starts:
					exon_Start_List.append(int(item))
				for item in exon_Ends:
					exon_End_List.append(int(item))
				exon_index_start = indexsearch(exon_Start_List, average)
				exon_index_end = indexsearch(exon_End_List, average)
				if exon_index_start != exon_index_end:
					exon_Number += 1
			else: 
				Intergenic_Number += 1 
				'''next = abs(average - TSS_List[min(TSSindex,(len(TSS_List)-1))])
				last = abs(average - TSS_List[max(0,(TSSindex-1))])
				if min(next, last) < 5000:
					if next < last: 
					else :
				print min(abs(average - TSS_List[min(TSSindex,(len(TSS_List)-1))]), abs(average - TSS_List[max(0,(TSSindex-1))])) < TSSthreshold'''
				#print chrom, islands[island_index].start, islands[island_index].end, islands[island_index].value
		Result['total'] = len(islands)
		Result['TSS'] = TSS_Number
		Result['GeneBody'] = GeneBody_Number
		Result['exon'] = exon_Number
		Result['Intergenic'] = Intergenic_Number
	return Result
		##print chrom, 'Total Islands:', len(islands), '. In TSS:', TSS_Number, '; In gene body region:', GeneBody_Number, '; In intergenic region:', Intergenic_Number




def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--bedfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="file with islands info in BED format")
	parser.add_option("-g", "--bedfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="file with summary info in bed graph format")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-w", "--window_size", action="store", type="int", dest="windowsize", help="window size of summary, negative if no summary", metavar="<int>")
	parser.add_option("-t", "--TSS_threshold", action="store", type="int", dest="TSSthreshold", help="upstream and downstream threshold of TSS region", metavar="<int>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	
	mm_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX', 'chrY', 'chrM']
	hg_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
	
	if opt.species == "mm8": 
		chroms = mm_chroms
	elif opt.species == "mm9": 
		chroms = mm_chroms
	elif opt.species == "hg18": 
		chroms = hg_chroms
	else:
		print 'Species Error'

	bed_vals = BED.BED(opt.species, opt.bedfile1, "BED3", 0)
	if opt.windowsize > 0:
		Histogram_bed = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH")
	else: Histogram_bed = 0
	
	coords = UCSC.KnownGenes(opt.known_genes)
	total_islands = 0
	total_TSS = 0
	total_GeneBody = 0
	total_exon = 0
	total_Intergenic = 0
	total_exon_length = 0.0
	total_gene_length = 0.0
	for chrom in chroms: 
		Result = get_island_regions(bed_vals, Histogram_bed, coords, chrom, opt.TSSthreshold, opt.windowsize)
		total_islands = total_islands + Result['total']
		total_TSS = total_TSS + Result['TSS']
		total_GeneBody = total_GeneBody + Result['GeneBody']
		total_exon = total_exon + Result['exon']
		total_Intergenic = total_Intergenic + Result['Intergenic']
		total_exon_length = total_exon_length + Result['exon_length']
		total_gene_length = total_gene_length + Result['gene_length']
	
	print 'Total Islands:', total_islands, '. In TSS:', total_TSS, '; In gene body region:', total_GeneBody, '(', total_exon, 'in exons and', total_GeneBody - total_exon, 'in introns); In intergenic region:', total_Intergenic, '; Gene length:', total_gene_length, '; Exon length:', total_exon_length
	#'(', total_exon_length / total_gene_length, ')', float(total_exon)/float(total_GeneBody)



if __name__ == "__main__":
	main(sys.argv)
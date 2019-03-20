#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
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
# wpeng@gwu.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from numpy import *

import BED
import GenomeData
import get_total_tag_counts
#import find_islands
#import Background_simulation_pr
import Background_island_probscore_statistics


def is_list_sorted(List):
        """
        Check if sorted in ascending order.
        input is a list of pure numbers.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(List)-1):
                if List[index] > List[index + 1]:
                        sorted = 0;
        return sorted;

def graining(List, win, step, score): 
	'''
	1 step coarse graining:
	List must be sorted!
	List (list) contains (start) coordinates of positive signals;
	min_win (int) is the window size in list, coarse graining will start from this resolution;
	step (int) is the number of windows in one graining unit;
	score (int) is the minimum number of positive elements in the graining unit to call the unit positive; 
	output is a library: "Count": list of positive unit number in each graining step;
	'''
	result = []
	endlimit = List[-1]
	i = 0
	k = 0
	while i <= endlimit and k < len(List):
		j = i + step * win
		h = k
		while h <= (len(List)-1) and List[h] < j:
			h += 1
		n = h - k
		if n >= score:
			result.append(i)
		k = h
		i = j
	return(result)


def coarsegraining(List, win_min, step, score):
	if (is_list_sorted(List) != 1):
		List.sort()
	Length_list = []
	Length_list.append(len(List))
	result_list = []
	result_list.append(List)
	win = win_min
	while len(List) > 0:
		print len(List)
		List = graining(List, win, step, score)
		Length_list.append(len(List))
		if len(List) > 0:
			result_list.append(List)
		win = win * step
	return Length_list, result_list


def union_islands_to_list(islandlist):
	'''input islandlist and output list are both lists of BED island objects'''
	islandlist.sort(key=operator.attrgetter('start'));
	List = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			List.append(current)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	List.append(current)
	return List

def write_islandlist(List, win):
	'''input a start list and universal island width, output a islandlist of BED objects 
	object.start = List[i]
	object.end = List[i] + win - 1'''
	output_list = []
	for item in List:
		output_list.append(BED.BED3('', item, item + win - 1))
	output_list.sort(key=operator.attrgetter('start'))
	return output_list


def backstep(islandlist, List, win):
	'''one step trace back'''
	#result_list = []
	#fine_islands = []
	addtional_islands = write_islandlist(List, win)
	'''for item in islandlist:
		start_left = (item.start - win) in List
		start_right = item.start in List
		if start_left and start_right:
			item.start = item.start - win
		elif (not start_left) and (not start_right):
			item.start = item.start + win
		end_left = (item.end + 1 - win) in List
		end_right = (item.end + 1) in List
		if end_left and end_right:
			item.end = item.end + win
		elif (not end_left) and (not end_right):
			item.end = item.end - win
		assert item.start < item.end'''
	return union_islands_to_list(islandlist + addtional_islands)
	
	
	

def stepback(start_list, end_list, List, win):
	'''
	one step trace back, refine islands defined by start_list and end_list into resolution win using List (of resolution win), plus the vanished islands in List
	win is the window size in List
	'''
	assert len(start_list) == len(end_list)
	result_start_list = []
	result_end_list = []
	i = 0
	j = 0
	while i < len(List) and j < len(start_list):
		if List[i] + win < start_list[j]:
			result_start_list.append(List[i])
			result_end_list.append(List[i] + win - 1)
			i += 1
		else:
			if List[i] + win == start_list[j]: # start left = 1
				if List[i+1] == start_list[j]: # start right = 1
					start = List[i] # 11 -> expand
				else: 
					start = start_list[j] # 10 -> keep
			elif  List[i] == start_list[j]: # start right = 1
				start = start_list[j] # 01 -> keep
			else:
				start = start_list[j] + win # 00 -> chrink
			while List[i] < end_list[j]:
				i += 1
			if List[i] == end_list[j] + 1: #right=1
				if List[i] == List[i-1] + win: #left=1
					end = List[i] + win - 1 # expand
				else:
					end = end_list[j] # 01 -> keep
			elif List[i-1] + win == end_list[j] + 1: #left = 1, 
				end = end_list[j] # 10 -> keep
			else:
				end = end_list[j] - win # 00 -> shrink
			result_start_list.append(start)
			result_end_list.append(end)
			j += 1
			i += 1
	return result_start_list, result_end_list
		


def traceback(List, win_min, step, level=0):
	'''
	Input is a list of lists. 
	'''
	win = win_min * pow(step, len(List)-1)
	islandlist = write_islandlist(List[-1], win)
	#start_list = List[-1]
	#end_list = []
	#for i in range(0, len(start_list)):
	#	end_list.append(start_list[i] + win - 1)
	i = 1
	while i < len(List)-level: 
		print len(islandlist)
		backlist = List[-i-1]
		win = win/step
		#starts = start_list
		#ends = end_list
		#(start_list, end_list) = stepback(starts, ends, backlist, win)
		islands = islandlist
		islandlist = backstep(islands, backlist, win)
		i += 1
	return islandlist

	
def main(argv):
	'''
	Coarse graining test chr1, input must only have chr1
	
	'''
	parser = OptionParser()
	
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summarygraph", action="store",type="string", dest="summarygraph", help="summarygraph", metavar="<file>")
	parser.add_option("-w", "--window_size(bp)", action="store", type="int", dest="window_size", help="window_size(in bps)", metavar="<int>")
	parser.add_option("-g", "--graining_size", action="store", type="int",  dest="step", help="graining unit size (>0)", metavar="<int>")
	parser.add_option("-e", "--score", action="store", type="int", dest="score", help="graining criterion, 0<score<= step", metavar="<int>")
	parser.add_option("-t", "--mappable_faction_of_genome_siz ", action="store", type="float", dest="fraction", help="mapable fraction of genome size", metavar="<float>")
	parser.add_option("-f", "--output_file", action="store", type="string", dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)

	print "Coarse graining: chr1 test";
	if opt.species in  GenomeData.species_chroms.keys():
		print "Species: ", opt.species;
		print "Window_size: ", opt.window_size;
		print "Coarse graining step: ", opt.step;
		print "Coarse graining score:", opt.score;
		
		total_read_count = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summarygraph);
		print "Total read count:", total_read_count
		genome_length = GenomeData.species_chrom_lengths[opt.species]['chr1']
		genome_length = int(opt.fraction * genome_length);

		average = float(total_read_count) * opt.window_size/genome_length; 
		print "Genome Length: ", genome_length;
		print "window average:", average;
		
		window_pvalue = 0.20;
		bin_size = 0.001;
		print "Window pvalue:", window_pvalue;
		background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_read_count, opt.window_size, 3 * opt.window_size, window_pvalue, genome_length, bin_size);
		min_tags_in_window = background.min_tags_in_window
		print "Minimum num of tags in a qualified window: ", min_tags_in_window
		
		print "Generate preprocessed data list"; 
		#read in the summary graph file
		bed_val = BED.BED(opt.species, opt.summarygraph, "BED_GRAPH");
		#generate the probscore summary graph file, only care about enrichment
		chrom = 'chr1' 
		if len(bed_val[chrom]) > 0:
			eligible_start_list = []
			for index in xrange(len(bed_val[chrom])):
				read_count = bed_val[chrom][index].value;
				if read_count >= min_tags_in_window:
					eligible_start_list.append(bed_val[chrom][index].start)
			print "Coarse graining:";
			(result_list, island_list) = coarsegraining(eligible_start_list, opt.window_size, opt.step, opt.score)
			print "Comnbing, No traceback...", len(island_list)
			islands = traceback(island_list, opt.window_size, opt.step, 0)
			print len(islands),"islands found!"
			f = open(opt.out_file, 'w')
			f.write('track type=bedGraph name=' + opt.out_file + '\n')
			for i in range(0, len(islands)):
				f.write('chr1\t' + str(int(islands[i].start)) + '\t' + str(int(islands[i].end)) + '\t1\n')
			f.close()
		else: 
			print "input data error!"
	else:
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
import separate_by_chrom
from GenomeData import *

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();



def is_sorted(list):
	"""
	Check if sorted in ascending order.
	input is a list of BED with chrom, start, end and value.
	output: sorted =1 or 0
	"""
	sorted = 1; 
	for index in range(0, len(list)-1):
		if list[index].start > list[index+1].start:
			sorted = 0;
	return sorted; 


def window_on_island(island, window):
	"""
	decide if a window is on an island 
	"""
	if (window.start >= island.start and window.end <= island.end):
		return 1;
	else:
		return 0;


def find_isle_windows(islands_list, window_size):
	"""
	Use islands to get all the windows with all tag numbers being zero.
	islands 
	 
	"""
	window_list = [];
	for island in islands_list:
		start = island.start;
		end = start + window_size -1;
		while end <= island.end:
			window_list.append(BED.BED_GRAPH(island.chrom, start, end, 0));
			start += window_size;
			end = start + window_size -1;
	return window_list;
	

	
			
def filter_out_uncovered_windows(islands_list, all_windows_list, window_size):
	"""
	This module uses islands list to filter the windows list, so that only those windows on islands are recorded 
	
	Both input lists are lists of BED_GRAPH objects, each of BED GRAPH object has tributes of start, end, and value. 
	The islands MUST BE MADE out of the corresponding summary graph file !!
	
	Windows that are in the gaps are also counted in.
	
	Both lists need to be sorted, which comes natually out of graphs and islands.   They are with commensuate boundaries, as the islands are built out of the windows. 
	
	"""
	assert( is_sorted(islands_list) == 1 );
	assert( is_sorted(all_windows_list) == 1 ); 
	
	isle_windows_list = find_isle_windows(islands_list, window_size);
	isle_windows_start_list=[];
	for window in isle_windows_list:
		isle_windows_start_list.append(window.start)
		
	all_windows_start_list = [];
	for window in all_windows_list:
		all_windows_start_list.append(window.start)
	
	# start position should be in the all_windows_start_list
	for index in xrange(len(isle_windows_start_list)):
		item = isle_windows_start_list[index];
		position = bisect.bisect_left(all_windows_start_list, item);
		if bisect.bisect_right(all_windows_start_list, item) - position == 1:
			isle_windows_list[index].value = all_windows_list[position].value;
	return isle_windows_list


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-b", "--islandfile", action="store", type="string", dest="islandbedfile", metavar="<file>", help="island file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="filtered summary graph file")
	parser.add_option("-w", "--window_size", action="store", type="int",  dest="window_size", help="window size of summary graph", metavar="<int>")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	summary_graph_extension=".summarygraph"
	island_extension=".islands"
	separate_by_chrom.separateByChrom(chroms, opt.bedfile, summary_graph_extension)
	separate_by_chrom.separateByChrom(chroms, opt.islandbedfile, island_extension)
	
	f = open(opt.out_file, 'w')
	for chrom in chroms:
		if separate_by_chrom.fileExists(chrom+summary_graph_extension) and separate_by_chrom.fileExists(chrom+island_extension):
			summary = BED.BED(opt.species, chrom+summary_graph_extension, "BED_GRAPH", 0);
			islands = BED.BED(opt.species, chrom+island_extension, "BED3", 0);
			result = filter_out_uncovered_windows(islands[chrom], summary[chrom], opt.window_size)
			for item in result:
				f.write(item.chrom + '\t' + str(item.start) +'\t'+ str(item.end) +'\t'+ str(item.value) + '\n')
	f.close()

	separate_by_chrom.cleanup(chroms, summary_graph_extension);
	separate_by_chrom.cleanup(chroms, island_extension);

if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
from GenomeData import *

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def filter_uncovered_windows(islands_list, windows_list):
	"""
	This module uses islands list to filter the windows list, so that only those windows on islands are recorded 
	
	Both input lists are lists of BED_GRAPH objects, each of BED GRAPH object has tributes of start, end, and value. 
	The islands MUST BE MADE out of the corresponding summary graph file !!
	
	Windows that are in the gaps are also counted in.
	
	Both lists are sorted.  They are with commensuate boundaries, as the islands are built out of the windows. 
	
	"""
	
	new_list=[];
	on_island_status = 0;
	island_index = 0;
	window_index = 0;
	
	while (island_index < len(islands_list) and window_index<len(windows_list)):
		island = islands_list[island_index];
		window = windows_list[window_index];
		
		if (window_on_island(island, window)==1):
				new_list.append(window);
				on_island_status = 1;
				window_index += 1;
		else:
			if (on_island_status == 1):
				island_index += 1;
				on_island_status = 0;
			else:
				window_index += 1;
	return new_list;


def window_on_island(island, window):
	"""
	decide if a window is on an island 
	"""
	if (window.start >= island.start and window.end <= island.end):
		return 1;
	else:
		return 0;


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-b", "--islandfile", action="store", type="string", dest="islandbedfile", metavar="<file>", help="island file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="filtered summary graph file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	summary = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0);
	islands = BED.BED(opt.species, opt.islandbedfile, "BED_GRAPH", 0);
	
	f = open(opt.out_file, 'w')
	for chrom in chroms:
		if chrom in summary.keys() and chrom in islands.keys():
			#print chrom
			result = filter_uncovered_windows(islands[chrom], summary[chrom])
			for item in result:
				f.write(chrom + '\t' + str(item.start) +'\t'+ str(item.end) +'\t'+ str(item.value) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)

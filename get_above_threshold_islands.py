#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED
import get_total_tag_counts
from species_chroms import *
from my_find_islands import *

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
			dest="species", help="mm8, hg18, background", metavar="<str>")
	parser.add_option("-i", "--islandfile", action="store", type="string",
			dest="islands_file",  metavar="<file>", 
			help="file with summary info in bed graph format")
	parser.add_option("-o", "--outfile", action="store", type="string",
			dest="outfile", help="output island file name", metavar="<file>")
	parser.add_option("-t", "--min_tag_count_in_island", action="store", type="int",
			dest="islands_minimum_tags", help="minimum number of tags required in each island", metavar="<int>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	if species_chroms.has_key(opt.species):
		total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.islands_file);
		print opt.islands_file;
		print "Total tag counts on original islands: "+ str(total_tag_counts);
		print "Island threshold is " + str(opt.islands_minimum_tags);
		
		find_region_above_threshold_from_file(opt.islands_file, opt.species, opt.islands_minimum_tags, opt.outfile);		
	else:
		print opt.species + " is not in the species list";
	

if __name__ == "__main__":
    	main(sys.argv)

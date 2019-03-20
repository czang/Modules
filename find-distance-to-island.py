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
from species_chroms import *
import get_total_tag_counts
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

def find_closest_island_around_promoter(island_start_list, island_end_list, position):
	start_index = bisect.bisect_left(island_start_list, position)
	end_index = bisect.bisect_left(island_end_list, position)
	if start_index - end_index == 1:  # TSS in an island
		return 0
	else: 
		assert end_index == start_index
		index = end_index
		if index == 0: 
			return abs(island_start_list[index] - position)
		elif index == len(island_end_list):
			return abs(position - island_end_list[index-1])
		else:
			up_distance = abs(position - island_end_list[index-1])
			down_distance = abs(island_start_list[index] - position)
			return min(up_distance, down_distance)


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island bed file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--known_genes_file", action="store", type="string", dest="bedfile2", metavar="<file>", help="file with candidate regions in BED file, calculating the center")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED3", 0); 
	candidate = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0)
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	#coords = UCSC.KnownGenes(opt.known_genes);
	#occupied_genes = {}; 

	
	for chrom in chroms:
		if chrom in candidate.keys():
			if chrom in bed_vals.keys():
				if len(bed_vals[chrom]) != 0:
					island_start_list=[];
					island_end_list=[];
					for item in bed_vals[chrom]:
						island_start_list.append(item.start);
						island_end_list.append(item.end);
					for item in candidate[chrom]:
						position = (item.start + item.end)/2
						item.value = find_closest_island_around_promoter(island_start_list, island_end_list, position)
			else:
				for item in candidate[chrom]:
					item.value = "9999999999"
	f = open(opt.out_file, 'w')
	for chrom in chroms:
		if chrom in candidate.keys():
			for item in candidate[chrom]:
				f.write(chrom + '\t' + str(item.start) + '\t' + str(item.end) + '\t' + str(item.value) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
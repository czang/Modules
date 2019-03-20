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

def find_closest_island_around_position(island_start_list, island_end_list, position):
	start_index = bisect.bisect_left(island_start_list, position)
	end_index = bisect.bisect_right(island_end_list, position)
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


def list2pdf(List, bin):
	histogram = [0]
	for score in List:
		index = int(score/bin);
		if index >= len(histogram):
			histogram += [0]*(index-len(histogram)+1);
		histogram[index] += 1.0;
	assert sum(histogram) == len(List)
	return histogram


def pdf2cdf(List):
	cdf = []
	n = 0
	cdf.append(n)
	for item in List:
		n += item
		cdf.append(n)
	return cdf


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="target island bed file")
	parser.add_option("-b", "--island_file", action="store", type="string", dest="bedfile2", metavar="<file>", help="file with candidate regions in BED file, calculating both boundaries")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for island and boundary distances")
	parser.add_option("-d", "--cdffile", action="store", type="string", dest="cdf_file", metavar="<file>", help="output file name for distance cdf")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED3", 0); 
	candidate = BED.BED(opt.species, opt.bedfile2, "BED3", 0)
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	#coords = UCSC.KnownGenes(opt.known_genes);
	#occupied_genes = {}; 

	f = open(opt.out_file, 'w')
	List = []
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
						left_dis = find_closest_island_around_position(island_start_list, island_end_list, item.start)
						right_dis = find_closest_island_around_position(island_start_list, island_end_list, item.end)
						List.append(left_dis)
						List.append(right_dis)
						f.write(chrom + '\t' + str(item.start) + '\t' + str(item.end) + '\t' + str(left_dis) + '\t' + str(right_dis) + '\n')
			else:
				for item in candidate[chrom]:
					List.append(9999999999)
					List.append(9999999999)
					f.write(chrom + '\t' + str(item.start) + '\t' + str(item.end) + '\t' + '9999999999' + '\t' + '9999999999' + '\n')
	f.close()
	
	print len(List)
	bin = 1000
	pdf = list2pdf(List, bin)
	cdf = pdf2cdf(pdf)
	f = open(opt.cdf_file, 'w')
	x = 0
	for item in cdf:
		f.write(str(x) + '\t' + str(item) + '\n')
		x += bin
	f.close()
	


if __name__ == "__main__":
	main(sys.argv)
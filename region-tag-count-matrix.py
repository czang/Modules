#!/usr/bin/env python
# Copyright (c) 2008 NHLBI, NIH
# Authors: Dustin Schones, Chongzhi Zang and Keji Zhao
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
# schonesde@mail.nih.gov).


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import bisect

import BED
import UCSC
import GenomeData
import SeparateByChrom
import associate_tags_with_regions


"""
The plan here is to get a profile for each modification across sets of
genes.

"""

upstream_dist = 25;
downstream_dist = 25;
intergenic_window = 5;
genic_percent_window = 0.05;
window_size = 2

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


chrom_length = {};
def getChromLengths(sp):
	for i in GenomeData.species_chrom_lengths[sp].keys():
		chrom_length[i] = GenomeData.species_chrom_lengths[sp][i]


def countTagsInWindow(start, end, tag_starts):
	start_index = bisect.bisect_left(tag_starts, start)
	end_index = bisect.bisect_left(tag_starts, end)
	return end_index - start_index


def getGenicCounts_fix_window(gene_coords, bed_file, shift, window_size, species):
	result = {}
	num_steps = float(downstream_dist) / float(intergenic_window)
	chroms = gene_coords.keys()
	SeparateByChrom.separateByChrom(chroms, bed_file, '.bed1')
	for chrom in gene_coords.keys():
		bed_vals = BED.BED(species, chrom+'.bed1', "BED2")
		plus_tag = []
		minus_tag = []
		for item in bed_vals[chrom]:
			if plus.match(item.strand):
				plus_tag.append(item.start + shift)
			elif minus.match(item.strand):
				minus_tag.append(item.start - 1 - shift)
		for g in gene_coords[chrom]:
			tag_count_list = []
			can_start_index = bisect.bisect_left(plus_tag, g.txStart - 100)
			can_end_index = bisect.bisect_right(plus_tag, g.txEnd + 100)
			can_tags_plus = plus_tag[can_start_index:(can_end_index + 1)]
			can_start_index = bisect.bisect_left(minus_tag, g.txStart - 100)
			can_end_index = bisect.bisect_right(minus_tag, g.txEnd + 100)
			can_tags_minus = minus_tag[can_start_index:(can_end_index + 1)]
			region_start_list = []
			region_end_list = []
			if plus.match(g.strand):
				start = g.txStart - upstream_dist    # upstream region
				end = start + intergenic_window
				region_start_list.append(start)
				region_end_list.append(end)
				for k in range(int(num_steps-1)):
					start += intergenic_window
					end += intergenic_window
					region_start_list.append(start)
					region_end_list.append(end)
				start = g.txStart    # genic region
				end = start + window_size
				region_start_list.append(start)
				region_end_list.append(end)
				while end < g.txEnd:
					start += window_size;
					end += window_size;
					region_start_list.append(start)
					region_end_list.append(end)
				start = g.txEnd    # downstream region
				end = start + intergenic_window
				region_start_list.append(start)
				region_end_list.append(end)
				for k in range(int(num_steps-1)):
					start += intergenic_window
					end += intergenic_window
					region_start_list.append(start)
					region_end_list.append(end)
			elif minus.match(g.strand):
				end = g.txEnd + upstream_dist    # upstream region
				start = end - intergenic_window
				region_start_list.append(start)
				region_end_list.append(end)
				for k in range(int(num_steps-1)):
					end -= intergenic_window
					start -= intergenic_window
					region_start_list.append(start)
					region_end_list.append(end)
				end = g.txEnd;    # genic region
				start = end - window_size;
				region_start_list.append(start)
				region_end_list.append(end)
				while start > g.txStart:
					start -= window_size;
					end -= window_size;
					region_start_list.append(start)
					region_end_list.append(end)
				start = g.txStart - intergenic_window    # downstream region
				end = start + intergenic_window
				region_start_list.append(start)
				region_end_list.append(end)
				for k in range(int(num_steps-1)):
					start -= intergenic_window
					end -= intergenic_window
					region_start_list.append(start)
					region_end_list.append(end)
			tag_plus_list = associate_tags_with_regions.find_readcount_on_regions(can_tags_plus, region_start_list, region_end_list)
			tag_minus_list = associate_tags_with_regions.find_readcount_on_regions(can_tags_minus, region_start_list, region_end_list)
			result[g.name] = tag_plus_list + tag_minus_list
	SeparateByChrom.cleanup(chroms, '.bed1')
	return result




def main(argv):
	parser = OptionParser()
	parser.add_option("-k", "--known_gene_file", action="store", type="string",
					  dest="genefile", help="file with known gene info", metavar="<file>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
					  dest="bedfile", help="file with tags in bed format", metavar="<file>")
	parser.add_option("-n", "--name", action="store", type="string",
					  dest="name", help="name for plotting", metavar="<str>")
	parser.add_option("-s", "--species", action="store", type="string",
					  dest="species", help="species", metavar="<str>")
	parser.add_option("-f", "--fragmentsize", action="store", type="int",
					  dest="fragmentsize", help="fragment size", metavar="<int>")


	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
		parser.print_help()
		sys.exit(1)

	gene_coords = UCSC.KnownGenes(opt.genefile)
	shift_length = int(round(opt.fragmentsize / 2))


	""" print everything out to a file """
	outfilename = '%s-scores' % opt.name
	outFile = open(outfilename, 'w')

	result = getGenicCounts_fix_window(gene_coords, opt.bedfile, shift_length, window_size, opt.species)
	for g in result.keys():
		List = result[g]
		strlist = []
		for i in range(0, len(List)):
			strlist.append(str(List[i]))
		outline = str(g) + '\t' + '\t'.join(strlist) + '\n'
		outFile.write(outline)
	outFile.close()


if __name__ == "__main__":
	main(sys.argv)


		

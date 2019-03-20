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


"""
The plan here is to get a profile for each modification across sets of
genes.

Break the 'gene region' up into:

-5kb ... -1 kb ...

...1st 5% of genes, 2nd 5% of gene,  .... last 5% of gene, ...

... +1 kb ... +5 kb

Get the tag densities in each region.
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

def getUpstreamCounts(gene_coords, bed_vals, shift):
	result = {}
	num_steps = float(upstream_dist) / float(intergenic_window)
	for chrom in bed_vals.keys():
		tag_starts = []
		for item in bed_vals[chrom]:
			if plus.match(item.strand):
				tag_starts.append(item.start + shift)
			elif minus.match(item.strand):
				tag_starts.append(item.start - 1 - shift)
		tag_starts.sort();

		if chrom in gene_coords.keys() and len(tag_starts)>0:
			for g in gene_coords[chrom]:
				tag_count_list = []
				if plus.match(g.strand):
					start = g.txStart - upstream_dist
					end = start + intergenic_window
					if start >= 0:
						tag_count_list.append(countTagsInWindow(start, end, tag_starts))
						for k in range(int(num_steps-1)):
							start += intergenic_window
							end += intergenic_window
							tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				elif minus.match(g.strand):
					end = g.txEnd + upstream_dist
					start = end - intergenic_window
					if end <= chrom_length[chrom]:
						tag_count_list.append(countTagsInWindow(start, end, tag_starts))
						for k in range(int(num_steps-1)):
							end -= intergenic_window
							start -= intergenic_window
							tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				result[g.name] = tag_count_list
	return result


def getDownstreamCounts(gene_coords, bed_vals, shift):
	result = {}
	num_steps = float(downstream_dist) / float(intergenic_window)
	for chrom in bed_vals.keys():
		tag_starts = []
		for item in bed_vals[chrom]:
			if plus.match(item.strand):
				tag_starts.append(item.start + shift)
			elif minus.match(item.strand):
				tag_starts.append(item.start - 1 - shift)
		tag_starts.sort();
		if chrom in gene_coords.keys() and len(tag_starts)>0:
			for g in gene_coords[chrom]:
				tag_count_list = []
				if plus.match(g.strand):
					start = g.txEnd
					end = start + intergenic_window
					if end <= chrom_length[chrom]:
						tag_count_list.append(countTagsInWindow(start, end, tag_starts))
						for k in range(int(num_steps-1)):
							start += intergenic_window
							end += intergenic_window
							tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				elif minus.match(g.strand):
					start = g.txStart - intergenic_window
					end = start + intergenic_window
					if start >= 0:
						tag_count_list.append(countTagsInWindow(start, end, tag_starts))
						for k in range(int(num_steps-1)):
							start -= intergenic_window
							end -= intergenic_window
							tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				result[g.name] = tag_count_list
	return result


def getGenicCounts_fix_window(gene_coords, bed_file, shift, window_size, species):
	result = {}
	num_steps = 1.0 / float(genic_percent_window)
	chroms = gene_coords.keys()
	SeparateByChrom.separateByChrom(chroms, bed_file, '.bed1')
	for chrom in gene_coords.keys():
		bed_vals = BED.BED(species, chrom+'.bed1', "BED2")
		tag_starts = []
		for item in bed_vals[chrom]:
			if plus.match(item.strand):
				tag_starts.append(item.start + shift)
			elif minus.match(item.strand):
				tag_starts.append(item.start - 1 - shift)
			tag_starts.sort()
		for g in gene_coords[chrom]:
			tag_count_list = []
			if plus.match(g.strand):
				start = g.txStart - upstream_dist    # upstream region
				end = start + intergenic_window
				tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				for k in range(int(num_steps-1)):
					start += intergenic_window
					end += intergenic_window
					tag_count_list.append(countTagsInWindow(start, end, tag_starts))
					print len(tag_count_list), tag_count_list[-1]
				start = g.txStart    # genic region
				end = start + window_size
				tag_count_list.append(countTagsInWindow(int(start), int(end), tag_starts))
				print len(tag_count_list), tag_count_list[-1]
				while end < g.txEnd:
					start += window_size;
					end += window_size;
					tag_count_list.append(countTagsInWindow(int(start), int(end), tag_starts))
					print len(tag_count_list), tag_count_list[-1]
				start = g.txEnd    # downstream region
				end = start + intergenic_window
				tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				print len(tag_count_list), tag_count_list[-1]
				for k in range(int(num_steps-1)):
					start += intergenic_window
					end += intergenic_window
					tag_count_list.append(countTagsInWindow(start, end, tag_starts))
					print len(tag_count_list), tag_count_list[-1]
			elif minus.match(g.strand):
				end = g.txEnd + upstream_dist    # upstream region
				start = end - intergenic_window
				tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				print len(tag_count_list), tag_count_list[-1]
				for k in range(int(num_steps-1)):
					end -= intergenic_window
					start -= intergenic_window
					tag_count_list.append(countTagsInWindow(start, end, tag_starts))
					print len(tag_count_list), tag_count_list[-1]
				end = g.txEnd;    # genic region
				start = end - window_size;
				tag_count_list.append(countTagsInWindow(int(start), int(end), tag_starts))
				print len(tag_count_list), tag_count_list[-1]
				while start > g.txStart:
					start -= window_size;
					end -= window_size;
					tag_count_list.append(countTagsInWindow(int(start), int(end), tag_starts))
					print len(tag_count_list), tag_count_list[-1]
				start = g.txStart - intergenic_window    # downstream region
				end = start + intergenic_window
				tag_count_list.append(countTagsInWindow(start, end, tag_starts))
				print len(tag_count_list), tag_count_list[-1]
				for k in range(int(num_steps-1)):
					start -= intergenic_window
					end -= intergenic_window
					tag_count_list.append(countTagsInWindow(start, end, tag_starts))
					print len(tag_count_list), tag_count_list[-1]
			result[g.name] = tag_count_list
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

	getChromLengths(opt.species);
	gene_coords = UCSC.KnownGenes(opt.genefile);
	#bed_vals = BED.BED(opt.species, opt.bedfile, "BED2")
	shift_length = int(round(opt.fragmentsize / 2))
	#num_tags = bed_vals.getNumVals();


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
	outFile.close();


if __name__ == "__main__":
	main(sys.argv)


		

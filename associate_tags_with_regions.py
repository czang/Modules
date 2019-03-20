#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
import GenomeData;
import utility
import separate_by_chrom

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

def tag_position(sline, fragment_size):
	if fragment_size >= 0:
		shift = int(round(fragment_size/2))
		if plus.match(sline[5]):
			return atoi(sline[1]) + shift
		elif minus.match(sline[5]):
			return atoi(sline[2]) - 1 - shift
	else:
		return int(round((atoi(sline[1]) + atoi(sline[2]))/2))

def find_readcount_on_regions(tag_position_list, region_start_list, region_end_list):
	'''The regions could overlap !!
	returns a list, with total tag number on this region, order as the region lists
	'''
	region_readcount_list = []
	assert len(region_start_list) == len(region_end_list)
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	
	for i in range(0, len(region_start_list)):
		if region_start_list[i] <= region_end_list[i]:
			s = bisect.bisect_left(tag_position_list, region_start_list[i])
			e = bisect.bisect_right(tag_position_list, region_end_list[i])
			tags = e-s;
			region_readcount_list.append(tags);
		else:
			region_readcount_list.append(0);
	return region_readcount_list


def find_read_sumscore_on_regions(tag_position_list, score_list, region_start_list, region_end_list):
	'''The regions could overlap !!
	tag position list and score list must be sorted!!
	'''
	region_readcount_list = []
	assert len(tag_position_list) == len(score_list)
	assert len(region_start_list) == len(region_end_list)
	if (utility.is_list_sorted(tag_position_list)==0):
		print "list not sorted, error will occur..."
	else: 
		for i in range(0, len(region_start_list)):
			if region_start_list[i] <= region_end_list[i]:
				s = bisect.bisect_left(tag_position_list, region_start_list[i])
				e = bisect.bisect_right(tag_position_list, region_end_list[i])
				score = 0.0
				for j in range(s,e):
					score += score_list[j]
				region_readcount_list.append(score)
			else:
				region_readcount_list.append(0.0)
	return region_readcount_list


def find_read_distance_weighted_score_on_regions(tag_position_list, site_list, d, f=4):
	'''The regions could overlap !!
	tag position list and score list must be sorted!!
	'''
	region_readcount_list = []
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	
	for i in range(0, len(site_list)):
		start = site_list[i] - d
		end = site_list[i] + d
		s = bisect.bisect_left(tag_position_list, start)
		e = bisect.bisect_right(tag_position_list, end)
		score = 0.0
		for j in range(s,e):
			x = abs(tag_position_list[j] - site_list[i])
			y = exp(-0.5 - f * (float(x)/100000))
			score += y
		region_readcount_list.append(score);
	return region_readcount_list


def find_strength_distance_weighted_score_on_regions(tag_position_list, score_list, site_list, d, f=4):
	'''The regions could overlap !!
	tag position list and score list must be sorted!!
	'''
	region_readcount_list = []
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	
	for i in range(0, len(site_list)):
		start = site_list[i] - d
		end = site_list[i] + d
		s = bisect.bisect_left(tag_position_list, start)
		e = bisect.bisect_right(tag_position_list, end)
		score = 0.0
		for j in range(s,e):
			x = abs(tag_position_list[j] - site_list[i])
			y = exp(-0.5 - f * (float(x)/100000)) * score_list[j]
			score += y
		region_readcount_list.append(score);
	return region_readcount_list


def rp(tag_position_list, site_list, start_list, end_list):
	'''The regions could overlap !!
	tag position list and score list must be sorted!!
	'''
	region_readcount_list = []
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	assert len(start_list) == len(end_list) == len(site_list)
	for i in range(0, len(start_list)):
		start = start_list[i]
		end = end_list[i]
		s = bisect.bisect_left(tag_position_list, start)
		e = bisect.bisect_right(tag_position_list, end)
		score = 0.0
		for j in range(s,e):
			x = abs(tag_position_list[j] - site_list[i])
			y = exp(-0.5 - 4 * (float(x)/100000))
			score += y
		region_readcount_list.append(score);
	return region_readcount_list


def rp_read_new(tag_position_list, site_list, start_list, end_list, lam = 10000):
	'''The regions could overlap !!
	tag position list and score list must be sorted!!
	'''
	region_readcount_list = []
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	assert len(start_list) == len(end_list) == len(site_list)
	for i in range(0, len(start_list)):
		start = start_list[i]
		end = end_list[i]
		s = bisect.bisect_left(tag_position_list, start)
		e = bisect.bisect_right(tag_position_list, end)
		score = 0.0
		for j in range(s,e):
			x = abs(tag_position_list[j] - site_list[i])
			expx = exp(-log(3.0) * float(x) / lam)
			y = 2 * expx / (1 + expx)
			score += y
		region_readcount_list.append(score);
	return region_readcount_list


def find_readdensity_on_regions(tag_position_list, region_start_list, region_end_list):
	'''The regions could overlap !!
	returns a list, with total tag number on this region, order as the region lists
	'''
	region_readcount_list = []
	assert len(region_start_list) == len(region_end_list)
	if (utility.is_list_sorted(tag_position_list)==0):
		tag_position_list.sort();
	
	for i in range(0, len(region_start_list)):
		if region_start_list[i] <= region_end_list[i]:
			s = bisect.bisect_left(tag_position_list, region_start_list[i])
			e = bisect.bisect_right(tag_position_list, region_end_list[i])
			tags = e-s;
			density = float(tags)/float((region_end_list[i] - region_start_list[i])*1000)
			region_readcount_list.append(density);
		else:
			region_readcount_list.append(0.0);
	return region_readcount_list

def find_readcount_on_islands(island_start_list, island_end_list, tag):
	"""
	Make sure the islands are sorted.
	
	islands are non-overlapping!
	"""
	
	index = bisect.bisect_right(island_start_list, tag);
	if index - bisect.bisect_left(island_end_list, tag) == 1:
		return index-1;
	else:
		return -1;
	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--rawreadfile", action="store", type="string", dest="readfile", metavar="<file>", help="raw read file in bed format")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="average size of a fragment after CHIP experiment")
	parser.add_option("-b", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>", help="island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	total = 0;
	
	islands = BED.BED(opt.species, opt.islandfile, "BED3", 0);
	if utility.fileExists(opt.readfile):
		separate_by_chrom.separateByChrom(chroms, opt.readfile, '.bed1');
	else:
		print opt.readfile, " not found";
		sys.exit(1)
	
	out = open(opt.out_file, 'w');
	for chrom in chroms:
		if chrom in islands.keys():
			island_list = islands[chrom];
			island_readcount_list=[0]*len(island_list);
			
			if utility.is_list_sorted(island_list) == 0:
				island_list.sort(key=operator.attrgetter('start'));
				
			island_start_list = []
			island_end_list = []
			for item in island_list:
				island_start_list.append(item.start)
				island_end_list.append(item.end)

			read_file = chrom + ".bed";
			f = open(read_file,'r')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					position = tag_position(sline, opt.fragment_size)
					index = find_readcount_on_islands(island_start_list, island_end_list, position);
					if index >= 0:
						island_readcount_list[index] += 1;
						total += 1;
			f.close();
							
			for index in xrange(len(island_list)):
				item = island_list[index];
				outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + str(island_readcount_list[index]) +"\n";	
				out.write(outline);		
							
	separate_by_chrom.cleanup(chroms, '.bed1');
	out.close();
	print "Total number of reads on islands are: ", total; 


if __name__ == "__main__":
	main(sys.argv)
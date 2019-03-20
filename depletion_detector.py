#!/usr/bin/env python
# Copyright (c) 2013 DFCI/HSPH
# Author: Chongzhi Zang
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
# czang@jimmy.harvard.edu)

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
from GenomeData import *
#import get_total_tag_counts
import SeparateByChrom
import associate_island_with_genes
import associate_tags_with_regions
import Utility

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def island_to_lists(island_start, island_end, window_size=200, step=20):
	start_list = []
	start = island_start
	end = start + window_size - 1
	while end <= island_end:
		start_list.append(start)
		start += step
		end = start + window_size - 1
	return start_list


def island_count_subset(island_start, island_end, position_list, count_list, step = 10):
	s = bisect.bisect_left(position_list, island_start)
	e = bisect.bisect_right(position_list, island_end)
	if s < e:
		positions = position_list[s:e]
		List = []
		start = positions[0]
		flag = position_list[s]
		while flag <= island_end:
			if flag in positions: 
				List.append(count_list[s])
				s +=1
			else: 
				List.append(0)
			flag += 10
		return List, start
	else:
		return [], -1


def score_list(start, count_list, w = 200, step = 10):
	l = len(count_list)
	avg = []
	for i in range(0,l-20):
		count = sum(count_list[i:i+20])
		avg.append(count)
	scores = []
	for i in range(0,l-60):
		score = (avg[i] + avg[i+40])/2 - avg[i+20]
		scores.append(max(score,0))
	return scores, start+300


def wig_to_lists(file):
	positionlist = []
	countlist = []
	f = open(file,'r')
	for line in f:
		if not (re.match("track", line) or re.match("variableStep",line)):
			line = line.strip()
			sline = line.split()
			positionlist.append(atoi(sline[0]))
			countlist.append(atoi(sline[1]))
	f.close()
	return positionlist, countlist


def peak_5(List):
	assert len(List) == 5
	if List[2] > List[1] and List[2] > List[3]:
		if List[1] > List[0] and List[3] > List[4]:
			return 1
		else:
			return -1
	else:
		return -1


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island file")
	parser.add_option("-g", "--wiggle", action="store", type="string", dest="wig", metavar="<file>", help="wig path")
	parser.add_option("-w", "--outwig", action="store", type="string", dest="out_wig", metavar="<file>", help="output wiggle file name")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output summit file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	if Utility.fileExists(opt.bedfile):
		SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1');
	else:
		print opt.bedfile, " not found";
		sys.exit(1)
	
	
	o = open(opt.out_file, 'w')
	for chrom in chroms:
		if Utility.fileExists(opt.wig+chrom+'.wig'):
			bed_vals = BED.BED(opt.species, chrom+'.bed1', "BED3", 0)
			islandlist = bed_vals[chrom]
			if len(islandlist) > 0:
				print 'reading', chrom, 'wig ...'
				(positionlist, countlist) = wig_to_lists(opt.wig+chrom+'.wig')
				print 'calculating...'
				xlist = []
				ylist = []
				for item in islandlist:
					(rawList, start) = island_count_subset(item.start, item.end, positionlist, countlist)
					if len(rawList) >= 60: 
						(scores, new_start) = score_list(start, rawList)
						for item in scores: 
							ylist.append(item)
						for i in range(0, len(scores)):
							xlist.append(new_start + i*10)
				assert len(xlist) == len(ylist)
				print 'writing wig track for', chrom
				f = open(chrom+'.out', 'w')
				f.write('track type=wiggle_0\n')
				f.write('variableStep chrom='+chrom+' span=10\n')
				for i in range(0,len(xlist)):
					if ylist[i] > 0:
						f.write(str(xlist[i]) + '\t' + str(ylist[i]) + '\n')
				f.close()
				print 'writing summits for', chrom
				for i in range(0,len(ylist)-5):
					if xlist[i+4] - xlist[i] < 50: #to confirm the 5 elements are consecutive 
						List = ylist[i:i+5]
						if peak_5(List) == 1:
							start = xlist[i+2] - 1
							end = start + 1
							score = List[2]
							o.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\n')
	o.close()
	
	SeparateByChrom.combineAllGraphFiles(chroms, '.out', opt.out_wig)
	SeparateByChrom.cleanup(chroms, '.bed1');
	SeparateByChrom.cleanup(chroms, '.out');


if __name__ == "__main__":
	main(sys.argv)

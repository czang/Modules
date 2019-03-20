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


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island file")
	parser.add_option("-g", "--wiggle", action="store", type="string", dest="wig", metavar="<file>", help="wig path")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
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
	
	
	for chrom in chroms:
		if Utility.fileExists(chrom+'.bed1'):
			print 'reading', chrom, 'wig ...'
			(positionlist, countlist) = wig_to_lists(opt.wig+chrom+'.wig')
			bed_vals = BED.BED(opt.species, chrom+'.bed1', "BED3", 0)
			islandlist = bed_vals[chrom]
			xlist = []
			ylist = []
			for item in islandlist:
				(rawList, start) = island_count_subset(item.start, item.end, positionlist, countlist)
				(scores, new_start) = score_list(start, rawList)
				for item in scores: 
					ylist.append(item)
				for i in range(0, len(scores)):
					xlist.append(new_start + i*10)
			print 'writing', chrom, 'out'
			f = open(chrom+'.out', 'w')
			f.write('track type=wiggle_0\n')
			f.write('variableStep chrom='+chrom+' span=10\n')
			assert len(xlist) == len(ylist)
			for i in range(0,len(xlist)):
				if ylist[i] > 0:
					f.write(str(xlist[i]) + '\t' + str(ylist[i]) + '\n')
			f.close()
	
	SeparateByChrom.combineAllGraphFiles(chroms, '.out', opt.out_file)
	SeparateByChrom.cleanup(chroms, '.bed1');
	#SeparateByChrom.cleanup(chroms, '.out');


if __name__ == "__main__":
	main(sys.argv)

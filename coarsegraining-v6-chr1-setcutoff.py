#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
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
# wpeng@gwu.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from numpy import *
import scipy.stats

import BED
import GenomeData
import get_total_tag_counts
import Background_island_probscore_statistics
import SeparateByChrom


'''version 4: 3-phase coarse graining, take the phase that has most 1 to next step. '''

def is_list_sorted(List):
        """
        Check if sorted in ascending order.
        input is a list of pure numbers.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(List)-1):
                if List[index] > List[index + 1]:
                        sorted = 0;
        return sorted;


def start_list_correlation_r(List, win, r, chrom_length):
	'''List must be sorted'''
	assert is_list_sorted(List) == 1
	x = List[0]%win
	sum = 0
	n = 0
	while x + r <= chrom_length:
		n += 1
		if x in List and (x + r) in List:
			sum += 1
		x += win
	if n > 0:
		return float(sum)/float(n)
	else:
		return 0.0


def start_list_correlation_function(List, win, chrom_length):
	xlist = []
	ylist = []
	file = open("cr_"+str(win)+".txt", 'w')
	for i in range(0, min(3, chrom_length/win)):
		r = i * win
		c = start_list_correlation_r(List, win, r, chrom_length)
		xlist.append(i)
		ylist.append(c)
		file.write(str(i)+'\t'+str(c)+'\n')
	file.close()
	return (xlist, ylist)


def correlation_length_fit(xlist, ylist):
	assert len(xlist) == len(ylist)
	loglist = []
	s = ylist[0] * ylist[0]
	for i in range(0, len(ylist)):
		loglist.append(log(max(ylist[i] - s, 0.000000000001)))
	(a,b,r,stderr,x) = scipy.stats.linregress(xlist,loglist)
	return -1.0/a

def graining(List, win, step, score): 
	'''
	1 step coarse graining, phase considered:
	List must be sorted!
	List (list) contains (start) coordinates of positive signals;
	win (int) is the window size in list, coarse graining will start from this resolution;
	step (int) is the number of windows in one graining unit;
	score (int) is the minimum number of positive elements in the graining unit to call the unit positive; 
	output is a list of positive unit number in each graining step;
	'''
	result = []
	endlimit = List[-1]
	for p in range(0, step):
		tmp_result = []
		i = List[0] - p * win
		k = 0
		while i <= endlimit and k < len(List):
			j = i + step * win
			h = k
			while h <= (len(List)-1) and List[h] < j:
				h += 1
			n = h - k
			if n >= score:
				tmp_result.append(i)
			k = h
			i = j
		if len(tmp_result) > len(result):
			result = tmp_result
	return(result)


def coarsegraining(List, win_min, step, score, genome_length):
	if (is_list_sorted(List) != 1):
		List.sort()
	Length_list = []
	Length_list.append(len(List))
	result_list = []
	result_list.append(List)
	win = win_min
	while len(List) > 0:
		#(xlist, ylist) = start_list_correlation_function(List, win, genome_length)
		print len(Length_list)-1, len(List)#, correlation_length_fit(xlist, ylist)
		List = graining(List, win, step, score)
		Length_list.append(len(List))
		if len(List) > 0:
			result_list.append(List)
		win = win * step
	return Length_list, result_list


def union_islands_to_list(islandlist, win):
	'''input islandlist and output list are both lists of BED island objects'''
	islandlist.sort(key=operator.attrgetter('start'));
	List = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end + 1 + win:
			List.append(current)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	List.append(current)
	return List


def write_islandlist(List, win):
	'''input a start list and universal island width, output a islandlist of BED objects 
	object.start = List[i]
	object.end = List[i] + win - 1'''
	output_list = []
	for item in List:
		output_list.append(BED.BED3('', item, item + win - 1))
	output_list.sort(key=operator.attrgetter('start'))
	return output_list


def backstep(islandlist, List, win):
	'''one step trace back'''
	#result_list = []
	#fine_islands = []
	addtional_islands = write_islandlist(List, win)
	for item in islandlist:
		start_left = (item.start - win) in List
		start_right = item.start in List
		if start_left and start_right:
			item.start = item.start - win
		elif (not start_left) and (not start_right):
			item.start = item.start + win
		end_left = (item.end + 1 - win) in List
		end_right = (item.end + 1) in List
		if end_left and end_right:
			item.end = item.end + win
		elif (not end_left) and (not end_right):
			item.end = item.end - win
		assert item.start < item.end
	return union_islands_to_list(islandlist + addtional_islands, win)


def traceback(List, win_min, step, level, genome_length):
	'''
	Input is a list of lists. 
	'''
	win = win_min * pow(step, len(List)-1)
	islandlist = write_islandlist(List[-1], win)
	backlist = List[-1]
	(xlist, ylist) = start_list_correlation_function(backlist, win, genome_length)
	correlation_length = correlation_length_fit(xlist, ylist)
	print len(backlist), correlation_length
	(xlist, ylist) = start_list_correlation_function(List[-2], win/step, genome_length)
	correlation_length_next = correlation_length_fit(xlist, ylist)
	print len(List[-2]), correlation_length_next
	i = 1
	while i < len(List)-level:
		backlist = List[-i-1]
		win = win/step
		if correlation_length > 1.0 and correlation_length_next > 1.0:
			print len(islandlist)
			islands = islandlist
			islandlist = backstep(islands, backlist, win)
		else:
			islandlist = write_islandlist(backlist, win)
			correlation_length = correlation_length_next
			(xlist, ylist) = start_list_correlation_function(List[-i-2], win/step, genome_length)
			correlation_length_next = correlation_length_fit(xlist, ylist)
			print len(List[-i-2]), correlation_length_next
		i += 1
	return islandlist


def main(argv):
	'''
	Coarse graining test chr1, input must only have chr1
	
	'''
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summarygraph", action="store",type="string", dest="summarygraph", help="summarygraph", metavar="<file>")
	parser.add_option("-w", "--window_size(bp)", action="store", type="int", dest="window_size", help="window_size(in bps)", metavar="<int>")
	parser.add_option("-g", "--graining_size", action="store", type="int",  dest="step", help="graining unit size (>0)", metavar="<int>")
	parser.add_option("-e", "--score", action="store", type="int", dest="score", help="graining criterion, 0<score<=graining_size", metavar="<int>")
	parser.add_option("-t", "--mappable_faction_of_genome_size", action="store", type="float", dest="fraction", help="mapable fraction of genome size", metavar="<float>")
	parser.add_option("-f", "--output_file", action="store", type="string", dest="out_file", help="output file name", metavar="<file>")
	parser.add_option("-m", "--min_tag", action="store", type="int",  dest="mintag", help="window minimum tag", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
        	parser.print_help()
        	sys.exit(1)

	print "Coarse-graining approach to identify ChIP-Seq enriched islands:"
	if opt.species in  GenomeData.species_chroms.keys():
		print "Species: ", opt.species;
		print "Window_size: ", opt.window_size;
		print "Coarse graining step: ", opt.step;
		print "Coarse graining score:", opt.score;
		chroms = ["chr1"]
		total_read_count = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summarygraph);
		print "Total read count:", total_read_count
		genome_length = GenomeData.species_chrom_lengths[opt.species]["chr1"];
		genome_length = int(opt.fraction * genome_length);

		average = float(total_read_count) * opt.window_size/genome_length; 
		print "Genome Length: ", genome_length;
		print "window average:", average;
		
		#window_pvalue = 0.20;
		bin_size = 0.001;
		#print "Window pvalue:", opt.window_pvalue;
		#background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_read_count, opt.window_size, 3 * opt.window_size, opt.window_pvalue, genome_length, bin_size);
		#min_tags_in_window = background.min_tags_in_window
		min_tags_in_window = opt.mintag
		#min_tags_in_window = 11
		print "Minimum num of tags in a qualified window: ", min_tags_in_window
		
		print "Generate preprocessed data list"; 
		#read in the summary graph file
		bed_val = BED.BED(opt.species, opt.summarygraph, "BED_GRAPH");
		#generate the probscore summary graph file, only care about enrichment
		for chrom in chroms: 
			if chrom in bed_val.keys() and len(bed_val[chrom]) > 0:
				eligible_start_list = []
				for index in xrange(len(bed_val[chrom])):
					read_count = bed_val[chrom][index].value;
					if read_count >= min_tags_in_window:
						eligible_start_list.append(bed_val[chrom][index].start)
				print "Coarse graining:";
				(result_list, island_list) = coarsegraining(eligible_start_list, opt.window_size, opt.step, opt.score, genome_length)
				print "Trace back...", len(island_list)
				islands = traceback(island_list, opt.window_size, opt.step, 0, genome_length)
				print len(islands), "islands found in", chrom
				f = open(chrom + ".islangstemp", 'w')
				for i in range(0, len(islands)):
					f.write(chrom + '\t' + str(int(islands[i].start)) + '\t' + str(int(islands[i].end)) + '\t1\n')
				f.close()
		o = open(opt.out_file, 'w')
		o.write('track type=bedGraph name=' + opt.out_file + '\n')
		o.close()
		SeparateByChrom.combineAllGraphFiles(chroms, ".islangstemp", opt.out_file)
		SeparateByChrom.cleanup(chroms, ".islangstemp")
		#else: 
			#print "input data error!"
	else:
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)
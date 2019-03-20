#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import separate_by_chrom

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def are_islands_sorted(island_list):
	"""
	Check if sorted in ascending order.
	input is a list of BED with chrom, start, end and value.
	output: sorted =1 or 0
	"""
	sorted = 1;
	for index in range(0, len(island_list)-1):
		if island_list[index].start> island_list[index+1].start:
			sorted = 0;
	return sorted;


def union_islands(islandlist):
	"""
	The islandlist MUST be pre-sorted according to start!!!
	"""
	start_list =[]
	end_list = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			start_list.append(current.start)
			end_list.append(current.end)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			i += 1
	start_list.append(current.start)
	end_list.append(current.end)
	return start_list, end_list


def region_overlap(start, end, start_list, end_list):
	"""
	Assuming non overlapping islands.
	"""
	assert (start <= end);
	start_position = bisect.bisect_right(end_list, start);
	end_position = bisect.bisect_left(start_list, end); 
	if (start_position < end_position): 
		return 1;
	else:
		return 0;


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--islandfile1", action="store", type="string", dest="islandfile1", metavar="<file>", help="file 1 with islands info to be compared")
	parser.add_option("-b", "--islandfile2", action="store", type="string", dest="islandfile2", metavar="<file>", help="file 2 with islands info to be compared")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-o", "--output", action="store", type="string", dest="overlapin1", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	total_overlap_number_1 = 0
	total_islands_1 = 0
	
	separate_by_chrom.separateByChrom(chroms, opt.islandfile1, '.island1')
	separate_by_chrom.separateByChrom(chroms, opt.islandfile2, '.island2')
	
	for chrom in chroms: 
		f = open(chrom + '.1in2', 'w')
		bed_vals_2 = BED.BED(opt.species, chrom+'.island2', "BED3", 0)
		if fileExists(chrom+'.island1') and len(bed_vals_2[chrom])>0:
			islandlist2 = bed_vals_2[chrom];
			if (are_islands_sorted(islandlist2) != 1):
				islandlist2.sort(key=operator.attrgetter('start'));
			(island2_start_list, island2_end_list) = union_islands(islandlist2)
			islands1 = open(chrom+'.island1', 'r')
			for line in islands1:
				if not re.match("#", line):
					total_islands_1 += 1
					line = line.strip()
					sline = line.split()
					start = int(sline[1])
					end = int(sline[2])
					name = (sline[3])
					if (region_overlap(start, end, island2_start_list, island2_end_list) == 1):
						f.write(name + '\t1\n')
						total_overlap_number_1 += 1;
					else:
						f.write(name + '\t0\n');
		elif fileExists(chrom+'.island1') and (len(bed_vals_2[chrom])==0):
			islands1 = open(chrom+'.island1', 'r')
			for line in islands1:
				if not re.match("#", line):
					total_islands_1 += 1
					line = line.strip()
					sline = line.split()
					f.write(sline[3] + '\t0\n');
		f.close()
	
	print "total number of island in "+opt.islandfile1+":     ", total_islands_1;
	print "total number of island in "+opt.overlapin1+":     ", total_overlap_number_1;
	
	separate_by_chrom.combineAllGraphFiles(chroms, '.1in2', opt.overlapin1);
	
	separate_by_chrom.cleanup(chroms, '.1in2')
	separate_by_chrom.cleanup(chroms, '.island1')
	separate_by_chrom.cleanup(chroms, '.island2')

if __name__ == "__main__":
	main(sys.argv)       
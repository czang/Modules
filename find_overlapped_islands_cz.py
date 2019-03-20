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


def write(item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	item is a BED3 object
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\n";
	out.write(outline);	


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--islandfile1", action="store", type="string", dest="islandfile1", metavar="<file>", help="file 1 with islands info to be compared")
	parser.add_option("-b", "--islandfile2", action="store", type="string", dest="islandfile2", metavar="<file>", help="file 2 with islands info to be compared")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-p", "--overlapin1", action="store", type="string", dest="overlapin1", metavar="<file>", help="file for islands in 1 overlapping with islands in 2")
	parser.add_option("-q", "--nonoverlapin1", action="store", type="string", dest="nonoverlapin1", help="file for islands in 1 not overlapping with islands in 2 ", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
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
		g = open(chrom + '.1notin2', 'w')
		bed_vals_1 = BED.BED(opt.species, chrom+'.island1', "BED3", 0)
		bed_vals_2 = BED.BED(opt.species, chrom+'.island2', "BED3", 0)
		if len(bed_vals_1[chrom]) > 0  and len(bed_vals_2[chrom])>0:
			islandlist1 = bed_vals_1[chrom];
			islandlist2 = bed_vals_2[chrom];
			total_islands_1 += len(bed_vals_1[chrom]);
			if (are_islands_sorted(islandlist1) != 1):
				islandlist1.sort(key=operator.attrgetter('start'));
			if (are_islands_sorted(islandlist2) != 1):
				islandlist2.sort(key=operator.attrgetter('start'));
			island2_start_list = []
			island2_end_list = [] 
			for item in islandlist2:
				island2_start_list.append(item.start);
				island2_end_list.append(item.end);
			for islandlist1_item in islandlist1:
				start = islandlist1_item.start;
				end = islandlist1_item.end;
				if (region_overlap(start, end, island2_start_list, island2_end_list) == 1):
					write (islandlist1_item, f);
					total_overlap_number_1 +=1;
				else:
					write (islandlist1_item, g);
		elif (len(bed_vals_1[chrom]) > 0) and (len(bed_vals_2[chrom])==0):
			total_islands_1 += len(bed_vals_1[chrom]);
			for islandlist1_item in bed_vals_1[chrom]:
				write (islandlist1_item, g);		
		f.close()
		g.close()
	
	print "total number of island in "+opt.islandfile1+":     ", total_islands_1;
	print "total number of island in "+opt.overlapin1+":     ", total_overlap_number_1;
	
	separate_by_chrom.combineAllGraphFiles(chroms, '.1in2', opt.overlapin1);
	separate_by_chrom.combineAllGraphFiles(chroms, '.1notin2', opt.nonoverlapin1);
	
	separate_by_chrom.cleanup(chroms, '.1in2')
	separate_by_chrom.cleanup(chroms, '.1notin2')
	separate_by_chrom.cleanup(chroms, '.island1')
	separate_by_chrom.cleanup(chroms, '.island2')

if __name__ == "__main__":
	main(sys.argv)       
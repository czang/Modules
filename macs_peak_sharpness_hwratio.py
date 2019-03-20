#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
import stats

import BED
from GenomeData import *
import separate_by_chrom


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to calculate the peak sharpness of a set of macs peaks. 
Peak sharpness is defined as summit height / average height of of the peak in the wiggle.

'''


def get_2lists_from_wig(wig_file):
	xlist = []
	ylist = []
	f = open(wig_file, 'r')
	for line in f:
		if (not re.match("track", line)) and (not re.match("variableStep", line)):
			line = line.strip()
			sline = line.split()
			assert len(sline) == 2
			xlist.append(int(sline[0]))
			ylist.append(float(sline[1]))
	f.close()
	return xlist, ylist


def calculate_sharpness(start, end, summit, wig_location_list, wig_count_list):
	'''standard deviation'''
	start_index = bisect.bisect_right(wig_location_list, start) - 1
	end_index = bisect.bisect_left(wig_location_list, end) 
	assert start_index < end_index
	if end_index - start_index > 1: 
		return stats.stdev(wig_count_list[start_index:end_index])
	else:
		return 0.0


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--peakfile", action="store", type="string", dest="bedfile", metavar="<file>", help="peak file with summits in Col4 in bedGraph format")
	parser.add_option("-w", "--wigfile", action="store", type="string", dest="wig", metavar="<file>", help="wigfile common part, chrom name will be auto added")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file in bedGraph format")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0)
	o = open(opt.out_file, 'w')
	for chrom in bed_vals.keys():
		(xlist, ylist) = get_2lists_from_wig(opt.wig+"_"+chrom+".wig")
		if len(bed_vals[chrom]) > 0:
			for island in bed_vals[chrom]:
				Summit = int(island.start + island.value - 1)
				value = calculate_sharpness(island.start, island.end, Summit, xlist, ylist)
				o.write(chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\t' + str(value) + '\n')
	o.close()
				

if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

## get BED module
import BED;
import UCSC;
from GenomeData import *;
import separate_by_chrom;


plus = re.compile("\+");
minus = re.compile("\-");


def get_regions(input_file):
	file = open(input_file, 'r')
	start_list = []
	end_list = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			start_list.append(atoi(sline[1]))
			end_list.append(atoi(sline[2]))
	file.close()
	start_list.sort()
	end_list.sort()
	return (start_list, end_list)


def get_list(input_file):
	file = open(input_file, 'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atoi(sline[1]))
	file.close()
	List.sort()
	return (start_list, end_list)


def separate_by_CTCF(start_list, end_list, CTCF_list):
	for CTCF in CTCF_list:
		index_start = bisect.bisect_left(start_list, CTCF)
		index_end = bisect.bisect_right(end_list, CTCF)
		if index_start != index_end:
			start_list.append(CTCF)
			end_list.append(CTCF)
			start_list.sort()
			end_list.sort()
	return (start_list, end_list)


def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-i", "--input_region_file", action="store", type="string", dest="region_file", metavar="<file>", help="input intergenic region file in BED format")
	parser.add_option("-c", "--input_CTCF_file", action="store", type="string", dest="CTCF_file", metavar="<file>", help="input CTCF binding sites file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
		chrom_lengths = species_chrom_lengths[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separate_by_chrom.separateByChrom(chroms, opt.region_file, '.bed1')
	separate_by_chrom.separateByChrom(chroms, opt.CTCF_file, '.bed2')
	Result = {}
	for chrom in chroms:
		region_filename = chrom + ".bed1"
		sites_filename = chrom + "bed2"
		if fileExists(region_filename):
			(start_list,end_list) = get_regions(region_filename)
			if fileExists(sites_filename):
				CTCF_list = get_list(sites_filename)
				Result[chrom] = separate_by_CTCF(start_list, end_list, CTCF_list)
			else:
				Result[chrom] = (start_list,end_list)
	separate_by_chrom.cleanup(chroms, '.bed1')
	separate_by_chrom.cleanup(chroms, '.bed2')
	file = open(opt.out_file, 'w')
	for chrom in Result.keys():
		(list1,list2) = Result[chrom]
		for i in range(0, len(list1)):
			file.write(chrom+'\t'+str(list1[i])+'\t'+str(list2[i])+'\n')
	file.close()


if __name__ == "__main__":
	main(sys.argv)
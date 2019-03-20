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

plus = re.compile("\+");
minus = re.compile("\-");


def find_gaps(list1, list2):
	assert len(list1) == len(list2)
	for i in range(1, len(list1)):
		assert list1[i-1] < list2[i-1]
		assert list2[i-1] < list1[i]
	start_list = list2
	end_list = list1
	del start_list[-1]
	del end_list[0]
	'''start_list.append(0)
	for item in list2:
		start_list.append(item)
	if length > list2[-1]:
		end_list.append(length)
	else:
		del start_list[-1]'''
	assert len(start_list)==len(end_list)
	checklist = []
	for i in range(0, len(start_list)):
		checklist.append(start_list[i])
		checklist.append(end_list[i])
	for i in range(1, len(checklist)):
		assert checklist[i] > checklist[i-1]
	return (start_list, end_list)


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-i", "--input_file", action="store", type="string", dest="in_file", metavar="<file>", help="input file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
		chrom_lengths = species_chrom_lengths[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals = BED.BED(opt.species, opt.in_file, "BED3", 0)
	#coords = UCSC.KnownGenes(opt.genes_file)
	Result = {}
	for chrom in chroms:
		if chrom in bed_vals.keys():
			island_list = bed_vals[chrom]
			island_list.sort(key=operator.attrgetter('start'))
			starts = []
			ends = []
			for item in island_list:
				starts.append(item.start)
				ends.append(item.end)
			#(region_starts, region_ends) = merge_overlap(gene_starts, gene_ends)
			Result[chrom] = find_gaps(starts, ends)
			#print chrom, len(coords[chrom]), len(region_starts)
	file = open(opt.out_file, 'w')
	for chrom in Result.keys():
		(list1,list2) = Result[chrom]
		for i in range(0, len(list1)):
			file.write(chrom+'\t'+str(list1[i])+'\t'+str(list2[i])+'\t'+str(list2[i]-list1[i]+1)+'\n')
	file.close()


if __name__ == "__main__":
	main(sys.argv)
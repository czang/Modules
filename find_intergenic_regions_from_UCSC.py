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


def merge_overlap(list1, list2):
	'''list1 must be sorted'''
	assert len(list1) == len(list2)
	start_list = []
	end_list = []
	start_index = 0
	end_index = 0
	while end_index < len(list2):
		assert start_index == end_index
		next_index = start_index + 1
		while next_index < len(list2) and list1[next_index] <= list2[end_index]:
			if list2[next_index] > list2[end_index]:
				end_index = next_index
			next_index += 1
		start_list.append(list1[start_index])
		end_list.append(list2[end_index])
		start_index = next_index
		end_index = next_index
	return (start_list, end_list)


def find_inter_regions(list1, list2, length):
	assert len(list1) == len(list2)
	for i in range(1, len(list1)):
		assert list1[i-1] < list2[i-1]
		assert list2[i-1] < list1[i]
	start_list = []
	end_list = list1
	start_list.append(0)
	for item in list2:
		start_list.append(item)
	if length > list2[-1]:
		end_list.append(length)
	else:
		del start_list[-1]
	assert len(start_list)==len(end_list)
	checklist = []
	for i in range(0, len(start_list)):
		checklist.append(start_list[i])
		checklist.append(end_list[i])
	for i in range(1, len(checklist)):
		assert checklist[i] > checklist[i-1]
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


def find_upstream_intergenic_region(gene_list, chrom, chrom_length, extension):
	result = []
	for i in range(0, len(gene_list)):
		g = gene_list[i]
		if plus.match(g.strand):
			start = g.txStart - extension
			if i == 0:
				result.append(g.name+'\t'+chrom+'\t-\t0\t'+str(start))
			elif start > (gene_list[i-1].txEnd + extension):
				result.append(g.name+'\t'+chrom+'\t-\t'+str(gene_list[i-1].txEnd + extension)+'\t'+str(start))
		elif minus.match(g.strand):
			start = g.txEnd + extension
			if i == len(gene_list)-1:
				result.append(g.name+'\t'+chrom+'\t+\t'+str(start)+'\t'+str(chrom_length))
			elif start < (gene_list[i+1].txStart-extension):
				result.append(g.name+'\t'+chrom+'\t+\t'+str(start)+'\t'+str(gene_list[i+1].txStart-extension))
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-i", "--input_genes_file", action="store", type="string", dest="genes_file", metavar="<file>", help="input genes file in UCSC known gene format")
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
	
	coords = UCSC.KnownGenes(opt.genes_file)
	Result = {}
	for chrom in chroms:
		if chrom in coords.keys():
			gene_list = coords[chrom]
			gene_list.sort(key=operator.attrgetter('txStart'))
			gene_starts = []
			gene_ends = []
			for item in gene_list:
				gene_starts.append(item.txStart)
				gene_ends.append(item.txEnd)
			(region_starts, region_ends) = merge_overlap(gene_starts, gene_ends)
			Result[chrom] = find_inter_regions(region_starts, region_ends, chrom_lengths[chrom])
			print chrom, len(coords[chrom]), len(region_starts)
	file = open(opt.out_file, 'w')
	for chrom in Result.keys():
		(list1,list2) = Result[chrom]
		for i in range(0, len(list1)):
			file.write(chrom+'\t'+str(list1[i])+'\t'+str(list2[i])+'\n')
	file.close()


if __name__ == "__main__":
	main(sys.argv)
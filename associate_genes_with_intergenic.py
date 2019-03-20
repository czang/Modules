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


def associate_genes_with_regions(start_list, end_list, genes, direction, limit):
	'''input data must be in the same chromosome, start_list and end_list must be sorted respectively'''
	Result = {}
	if minus.match(direction):
		for g in genes:
			if plus.match(g.strand):
				position = g.txStart
				index = bisect.bisect_left(end_list, position)
				if index > 0:
					index -= 1
					Result[g.name] = '-\t'+str(max(start_list[index], end_list[index] - limit))+'\t'+str(end_list[index])
			elif minus.match(g.strand):
				position = g.txEnd
				index = bisect.bisect_left(start_list, position)
				if index < len(start_list):
					Result[g.name] = '+\t'+str(start_list[index])+'\t'+str(min(end_list[index], start_list[index] + limit))
	elif plus.match(direction):
		for g in genes:
			if plus.match(g.strand):
				position = g.txEnd
				index = bisect.bisect_left(start_list, position)
				if index < len(start_list):
					Result[g.name] = '+\t'+str(start_list[index])+'\t'+str(min(end_list[index], start_list[index] + limit))
			elif minus.match(g.strand):
				position = g.txStart
				index = bisect.bisect_left(end_list, position)
				if index > 0:
					index -= 1
					Result[g.name] = '-\t'+str(max(start_list[index], end_list[index] - limit))+'\t'+str(end_list[index])
	else:
		print "Direction error, + or -"
	return Result


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
	parser.add_option("-g", "--input_gene_file", action="store", type="string", dest="gene_file", metavar="<file>", help="input gene file in UCSC format")
	parser.add_option("-d", "--direction", action="store", type="string", dest="direction", metavar="<str>", help="- for upstream, + for downstream")
	parser.add_option("-n", "--limit", action="store", type="int", dest="limit", metavar="<int>", help="limit for an intergenic region for a gene")
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
	
	coords = UCSC.KnownGenes(opt.gene_file);
	separate_by_chrom.separateByChrom(chroms, opt.region_file, '.bed')
	Result = {}
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom]
			region_filename = chrom + ".bed"
			if fileExists(region_filename):
				(start_list,end_list) = get_regions(region_filename)
				Result[chrom] = associate_genes_with_regions(start_list, end_list, genes, opt.direction, opt.limit)
	separate_by_chrom.cleanup(chroms, '.bed')
	file = open(opt.out_file, 'w')
	for chrom in Result.keys():
		regions = Result[chrom]
		for g in regions.keys():
			file.write(g+'\t'+chrom+'\t'+regions[g]+'\t0\t0\t0\t0\t0\n')
	file.close()


if __name__ == "__main__":
	main(sys.argv)
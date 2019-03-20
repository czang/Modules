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
import SeparateByChrom
import Utility

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def find_nearest_genes(summit, genes, max_distance):
	genes.sort(key=operator.attrgetter('txStart'))
	result_dic = {}
	plus_up = max_distance
	plus_down = max_distance
	minus_up = max_distance
	minus_down = max_distance
	for g in genes:
		if plus.match(g.strand):
			TSS = g.txStart
			distance = TSS - summit
			if distance >= 0:
				if distance < plus_up:
					plus_up_gene = g.name
					plus_up = distance
			else:
				distance = abs(distance)
				if distance < plus_down:
					plus_down_gene = g.name
					plus_down = distance
		elif minus.match(g.strand):
			TSS = g.txEnd
			distance = summit - TSS
			if distance >= 0:
				if distance < minus_up:
					minus_up_gene = g.name
					minus_up = distance
			else:
				distance = abs(distance)
				if distance < minus_down:
					minus_down_gene = g.name
					minus_down = distance
	if plus_up < max_distance:
		result_dic[plus_up_gene] = plus_up
	if plus_down < max_distance:
		result_dic[plus_down_gene] = -plus_down
	if minus_up < max_distance:
		result_dic[minus_up_gene] = minus_up
	if minus_down < max_distance:
		result_dic[minus_down_gene] = -minus_down
	return result_dic


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-d", "--max_distance", action="store", type="int", dest="max_distance", help="maximum distance to search for taget, default 100Kb", metavar="<int>", default=100000)
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species]
	else:
		print "This species is not recognized, exiting"
		sys.exit(1)
	
	if Utility.fileExists(opt.bedfile):
		SeparateByChrom.separateByChrom(chroms, opt.bedfile, '.bed1')
	else:
		print opt.bedfile, " not found"
		sys.exit(1)
	
	coords = UCSC.KnownGenes(opt.known_genes)
	all_result = {}
	for chrom in chroms:
		read_file = chrom + ".bed1"
		if Utility.fileExists(read_file) and chrom in coords.keys():
			f = open(read_file,'r')
			genes = coords[chrom]
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					summit = int(sline[2])
					result = find_nearest_genes(summit, genes, opt.max_distance)
					for g in result.keys():
						if g in all_result.keys():
							if abs(result[g]) < abs(all_result[g]):
								all_result[g] = result[g]
						else:
							all_result[g] = result[g]
			f.close()
	o = open(opt.out_file, 'w')
	for g in all_result.keys():
		o.write(g + '\t' + str(all_result[g]) + '\n')
	o.close()
	SeparateByChrom.cleanup(chroms, '.bed1')


if __name__ == "__main__":
	main(sys.argv)


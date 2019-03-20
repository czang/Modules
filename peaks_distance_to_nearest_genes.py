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


def find_nearest_genes(summit_file, TSS_file, max_distance=100000):
	TSS_dic = {}
	g = open(TSS_file, 'r')
	for line in g:
		line = line.strip()
		sline = line.split()
		if plus.match(sline[5]):
			TSS_dic[sline[3]] = atoi(sline[1])
		elif minus.match(sline[5]):
			TSS_dic[sline[3]] = atoi(sline[2])
	g.close()
	result_dic = {}
	f = open(summit_file, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			summit = int(sline[2])
			current_distance = max_distance
			for gene in TSS_dic.keys():
				distance = abs(summit - TSS_dic[gene])
				if distance < current_distance:
					current_distance = distance
					if gene in result_dic.keys():
						if distance < result_dic[gene]:
							result_dic[gene] = distance
					else:
						result_dic[gene] = distance
	f.close()
	return result_dic

def find_distance_to_nearest_genes(summit_file, TSS_file):
	TSS = []
	g = open(TSS_file, 'r')
	for line in g:
		line = line.strip()
		sline = line.split()
		if plus.match(sline[5]):
			TSS.append(atoi(sline[1]))
		elif minus.match(sline[5]):
			TSS.append(atoi(sline[2]))
	g.close()
	TSS.sort()
	result_dic = {}
	f = open(summit_file, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			summit = int(sline[2])
			tag = bisect.bisect_left(TSS, summit)
			if tag == 0:
				d = abs(TSS[tag] - summit)
			elif tag == len(TSS):
				d = abs(TSS[tag-1] - summit)
			else: 
				d = min(abs(TSS[tag] - summit), abs(TSS[tag-1] - summit))
			result_dic{sline[3]} = d
	f.close()
	return result_dic


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="ChIP seq read file")
	parser.add_option("-g", "--TSS_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with genes TSS in the format of name, chrom, TSS")
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
		print opt.bedfile, "not found"
		sys.exit(1)
	
	if Utility.fileExists(opt.known_genes):
		SeparateByChrom.separateByChrom(chroms, opt.known_genes, '.bed2')
	else:
		print opt.known_genes, "not found"
		sys.exit(1)
	

	all_result = {}
	for chrom in chroms:
		read_file = chrom + ".bed1"
		TSS_file = chrom + ".bed2"
		if Utility.fileExists(read_file) and Utility.fileExists(TSS_file):
			result = find_nearest_genes(read_file, TSS_file, opt.max_distance)
			for g in result.keys():
					all_result[g] = result[g]

	o = open(opt.out_file, 'w')
	for g in all_result.keys():
		o.write(g + '\t' + str(all_result[g]) + '\n')
	o.close()
	SeparateByChrom.cleanup(chroms, '.bed1')
	SeparateByChrom.cleanup(chroms, '.bed2')


if __name__ == "__main__":
	main(sys.argv)


#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from species_chroms import *
import get_total_tag_counts
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


'''This module is find 200bp windows that dont have any modification tags. '''


modifications = ['H2AK5ac','H2AK9ac','H2AK127me1','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K23me2','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K79me2','H4K91ac','H4R3me2']


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="path of summary graph files")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="windowsize", help="window size in bp", metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
		chrom_lengths = species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	result = {}
	for chrom in chroms:
		result[chrom] = []
		i = 0
		while i < chrom_lengths[chroms.index(chrom)]:
			result[chrom].append(i)
			i += opt.windowsize
	
	for mod in modifications:
		print mod
		bed_vals = {}
		bed_vals = BED.BED(opt.species, opt.bedfile+'/'+mod+'_summary.graph', "BED_GRAPH", 0);
		for chrom in chroms:
			if chrom in bed_vals.keys():
				for item in bed_vals[chrom]:
					if item.start in result[chrom]:
						result[chrom].remove(item.start)
	
	f = open(opt.out_file,'w')
	for chrom in chroms:
		for item in result[chrom]:
			f.write(chrom+'\t'+str(item)+'\t'+str(item+opt.windowsize-1)+'\t0\n')
	f.close()

if __name__ == "__main__":
	main(sys.argv)
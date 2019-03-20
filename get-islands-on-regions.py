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
import separate_by_chrom
import associate_island_with_genes


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="islandfile", metavar="<file>", help="file with islands info in BED format")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="summaryfile", metavar="<file>", help="file with summary graph info in BED format")
	parser.add_option("-k", "--genefile", action="store", type="string", dest="genefile", metavar="<file>", help="file with region info in UCSC known gene format")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="windowsize", metavar="<int>", help="window size in bps")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	separate_by_chrom.separateByChrom(chroms, opt.islandfile, '.island')
	separate_by_chrom.separateByChrom(chroms, opt.summaryfile, '.bed')
	coords = UCSC.KnownGenes(opt.genefile)
	
	result = {}
	for chrom in chroms:
		if chrom in coords.keys():
			if associate_island_with_genes.fileExists(chrom+'.island'):
				island_vals = BED.BED(opt.species, chrom+'.island', "BED_GRAPH", 0)
				summary_vals = BED.BED(opt.species, chrom+'.bed', "BED_GRAPH", 0)
				(island_position_list, island_value_list) = associate_island_with_genes.generate_island_position_and_value_lists(island_vals[chrom], summary_vals[chrom], opt.windowsize)
				region_names = []
				region_start_list = []
				region_end_list = []
				for g in coords[chrom]:
					region_names.append(g.name)
					region_start_list.append(g.txStart)
					region_end_list.append(g.txEnd)
				value_list = associate_island_with_genes.place_island_center_on_regions(island_position_list, island_value_list, region_start_list, region_end_list)
				result[chrom] = (region_names, value_list)
	
	f = open(opt.out_file, 'w')
	for chrom in result.keys():
		(List1, List2)= result[chrom]
		assert len(List1) == len(List2)
		for i in range(0, len(List1)):
			f.write(List1[i]+'\t'+str(List2[i])+'\n')
	f.close()
	
	separate_by_chrom.cleanup(chroms, '.island')
	separate_by_chrom.cleanup(chroms, '.bed')


if __name__ == "__main__":
	main(sys.argv)
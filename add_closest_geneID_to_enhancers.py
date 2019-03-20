#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
import UCSC
from gene_set_manipulation import *
from associate_island_with_genes import *
from species_chroms import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-i", "--infile", action="store", type="string", dest="bedfile", metavar="<file>", help="input bed summary graph file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", 0);
	for chrom in bed_vals.keys():
		if chrom not in chroms:
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes);
	
	f = open(opt.outfile, 'w')
	for chrom in bed_vals.keys():
		if chrom in coords.keys():
			genes = coords[chrom]
			TSS_list = get_TSS_list(genes)
			
			for item in bed_vals[chrom]:
				location = (item.start + item.end)/2
				index = associate_a_location_to_genes(location, TSS_list)
				f.write(genes[index].name+'\t'+chrom+'\t+\t'+str(item.start)+'\t'+str(item.end)+'\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
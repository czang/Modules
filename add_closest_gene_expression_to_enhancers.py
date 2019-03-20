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
import associate_binary_modification_with_expression

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
	
	enhancers = UCSC.KnownGenes(opt.bedfile)
	coords = UCSC.KnownGenes(opt.known_genes);
	expression = associate_binary_modification_with_expression.get_expression_data_dic(opt.known_genes, 12)
	
	f = open(opt.outfile, 'w')
	for chrom in enhancers.keys():
		if chrom in coords.keys():
			genes = coords[chrom]
			TSS_list = get_TSS_list(genes)
			
			for item in enhancers[chrom]:
				location = (item.txStart + item.txEnd)/2
				index = associate_a_location_to_genes(location, TSS_list)
				f.write(item.name+'\t'+chrom+'\t+\t'+str(item.txStart)+'\t'+str(item.txEnd)+'\t0\t0\t0\t0\t0\t0\t0\t'+str(expression[genes[index].name])+'\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
#!/usr/bin/env python
# Copyright (c) 2009 GWU & NHLBI, NIH
# Authors: Chongzhi Zang, Weiqun Peng and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# zang@gwmail.gwu.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from GenomeData import *
import get_total_tag_counts
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="island bed file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for genes and tag numbers")
	
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
	for chrom in chroms:
		if chrom not in bed_vals.keys(): 
			print chrom, " name is not the same as the stored one";
			sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes);
	occupied_genes = {}; 
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			bed_graph_list = bed_vals[chrom];
			current_genes = associate_island_with_genes.find_closest_island_around_promoter(bed_graph_list, genes);
			for gene in current_genes.keys():
				occupied_genes[gene] = current_genes[gene]
	
	f = open(opt.out_file, 'w')
	for g in occupied_genes.keys():
		f.write(g + '\t' + str(occupied_genes[g]) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
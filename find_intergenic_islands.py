#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Authors: Chongzhi Zang and Xiaole Shirley Liu
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
# czang@jimmy.harvard.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from GenomeData import *
from associate_island_with_genes import *

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="original islands file")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-p", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-o", "--outfileoccupiedgenes", action="store", type="string", dest="out_file", metavar="<file>", help="output file name for intergenic region islands")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	region_type = 'GenePromoter'
	
	bed_vals = BED.BED(opt.species, opt.bedfile, "BED3", 0); #islands 
	coords = UCSC.KnownGenes(opt.known_genes);
	total = 0;
	f = open(opt.out_file, 'w')
	#overlapped_islands = [];
	
	for chrom in chroms:
		if chrom in bed_vals.keys():
			bed_graph_list = bed_vals[chrom];
			if chrom in coords.keys():
				genes = coords[chrom];
				overlapped_islands = find_domains_not_on_genes(bed_graph_list, genes, region_type, opt.promoter_extension)
				for i in range(0, len(overlapped_islands)):
					if overlapped_islands[i] == 0:
						island = bed_graph_list[i]
						f.write(island.chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\n')
						total+=1;
			else:
				for island in bed_graph_list:
					f.write(island.chrom + '\t' + str(island.start) + '\t' + str(island.end) + '\n')
					total+=1;
	f.close()
	print total, "islands of", opt.bedfile, "found in intergenic region."


if __name__ == "__main__":
	main(sys.argv)

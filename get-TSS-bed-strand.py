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

import UCSC;
from GenomeData import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-t", "--promoter_threshold", action="store", type="int", dest="TSSthreshold", help="upstream and downstream threshold of promoter region", metavar="<int>")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);

	coords = UCSC.KnownGenes(opt.known_genes)
	f = open(opt.out_file, 'w')
	for chrom in chroms: 
		if chrom in coords.keys():
			genes = coords[chrom]
			for g in genes: 
				if plus.match(g.strand):
					TSS = g.txStart
				elif minus.match(g.strand):
					TSS = g.txEnd
				start = TSS - opt.TSSthreshold
				end = TSS + opt.TSSthreshold
				f.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + g.name + '\t0\t' + g.strand + '\n')
	f.close()
	


if __name__ == "__main__":
	main(sys.argv)
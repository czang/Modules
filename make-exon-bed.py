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

## get BED module
import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	coords = UCSC.KnownGenes(opt.known_genes)
	f = open(opt.out_file, 'w')
	for chrom in coords.keys():
		genes = coords[chrom]
		for gene in genes:
			exon_Starts = gene.exonStarts.replace(',',' ')
			exon_Ends = gene.exonEnds.replace(',',' ')
			exon_Starts = exon_Starts.split()
			exon_Ends = exon_Ends.split()
			assert len(exon_Starts) == len(exon_Ends) == int(gene.exonCount)
			for i in range(0, len(exon_Starts)):
				f.write(chrom + '\t' + exon_Starts[i] + '\t' + exon_Ends[i] + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
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


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="windowsize", help="window size in bps", metavar="<int>")
	parser.add_option("-p", "--windowpvalue", action="store", type="float", dest="pvalue", help="P value, only windows those have a probability of having such number of tags is less than p are counted", metavar="<float>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	genomelength = 0.0
	for item in species_chrom_lengths[opt.species]:
		genomelength += item
	total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile)
	
	f = open(opt.out_file, 'w')
	for i in range(1,15):
		average = total_tag_counts / genomelength * i * float(opt.windowsize)
		threshold = associate_island_with_genes.find_threshold(opt.pvalue, average)
		f.write(str(i)+'\t'+str(threshold)+'\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
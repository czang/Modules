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
from get_total_tag_counts import *
from associate_island_with_genes import *

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

''''''

modifications_38 = ['H2AK5ac','H2AK9ac','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K91ac']


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="path of summary graph files")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-t", "--startextension", action="store", type="int", dest="start_extension", help="upstream threshold of gene region", metavar="<int>")
	parser.add_option("-e", "--endextension", action="store", type="int", dest="end_extension", help="downstream threshold of gene region", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.known_genes)
	histogram = [0]*10000
	total_data = 0
	
	for mod in modifications_38:
		print mod
		bed_vals = {}
		bed_vals = BED.BED(opt.species, opt.bedfile+'/'+mod+'_summary.graph', "BED_GRAPH", 0)
		total = get_total_tag_counts_bed_graph(opt.bedfile+'/'+mod+'_summary.graph')
		for chrom in chroms:
			if chrom in bed_vals.keys() and chrom in coords.keys():
				print '\t',chrom
				List = place_summary_on_enhancers(bed_vals[chrom], coords[chrom], opt.start_extension, opt.end_extension, 200)
				for item in List:
					histogram[int(item/float(total)*10000000.0)] += 1
					total_data += 1
	
	f = open('enhancer_window_histogram_all_10M.txt','w')
	for i in range(0, len(histogram)):
		f.write(str(i) + '\t' + str(float(histogram[i])/float(total_data)*100.0) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
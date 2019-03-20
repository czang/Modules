#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import get_total_tag_counts
import expression_histogram
import list_histogram
from species_chroms import *

Dir = os.getcwd();


modifications_38 = ['H2AK5ac','H2AK9ac','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K91ac']


def get_all_data_list(summary_file):
	all_list = [0]*10000
	for mod in modifications_38:
		total = get_total_tag_counts.get_total_tag_counts_bed_graph(summary_file+'/'+mod+'_summary.graph')
		data_list = list_histogram.get_data_list(summary_file+'/'+mod+'_summary.graph',3)
		for item in data_list:
			all_list[int(item/float(total)*10000000.0)] += 1
	return all_list


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="path of summary files")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_all_data_list(opt.infile)
	
	genomelength = 0.0
	for item in species_chrom_lengths['hg18']:
		genomelength += item
	
	f = open(opt.output,'w')
	for i in range(0, len(List)):
		f.write(str(i) + '\t' + str(float(List[i])/(38*genomelength/200)*100) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
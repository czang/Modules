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
import gene_set_manipulation
from associate_binary_modification_with_expression import *


plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

'''This module is to calculate the modification data specially on gene promoters. Modification data are gained from summary graph file. Given a p-value, if total tag count on the promoter region is more than the threshold, the number is recorded with the gene. Otherwise record 0.'''

Backbone_17 = ['H2A_Z','H2BK5ac','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K18ac','H3K27ac','H3K36ac','H4K5ac','H4K8ac','H4K91ac']

Modifications = Backbone_17

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--summarygraphfile", action="store", type="string", dest="tagfile", metavar="<file>", help="summary graph file")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-t", "--total_tag_number_file", action="store", type="string", dest="total_tags", metavar="<file>", help="file with total tag numbers")
	#parser.add_option("-o", "--out_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	totals = get_expression_data_dic(opt.total_tags, 1)
	for mod in Modifications:
		assert mod in totals.keys()
	gene_list = gene_set_manipulation.get_gene_list(opt.known_genes, 0)
	result = {}
	for gene in gene_list:
		result[gene] = 0.0
	for mod in Modifications:
		tags = get_expression_data_dic(mod+opt.tagfile, 1)
		for gene in gene_list:
			assert gene in tags.keys()
			tag = tags[gene] / totals[mod] *1000000
			result[gene] += tag
	List = []
	for gene in result.keys():
		List.append(result[gene])
	print average(List)


if __name__ == "__main__":
	main(sys.argv)
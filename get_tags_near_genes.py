#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


import BED
import UCSC
from GenomeData import *
#import get_total_tag_counts
import SeparateByChrom
import associate_island_with_genes
import associate_tags_with_regions
import Utility

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-b", "--bedfilepath", action="store", type="string", dest="bedpath", metavar="<file>", help="ChIP seq read file path")
	parser.add_option("-g", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes in UCSC format")
	parser.add_option("-r", "--'Promoter', 'GeneBody', 'GenePromoter', or 'GeneExtended'", action="store", type="string", dest="region_type", metavar="<str>", help="region to count tags in", default="Promoter")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>", default=1000)
	parser.add_option("-o", "--outfilepath", action="store", type="string", dest="out_file_path", metavar="<file>", help="output files directory", default="./")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	allowed_region_type = ['Promoter', 'GeneBody', 'GenePromoter', 'GeneExtended'];
	if opt.region_type not in allowed_region_type: 
		print "The region type is not recognized, exiting";
		sys.exit(1);
	
	
	coords = UCSC.KnownGenes(opt.known_genes);
	
	
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom];
			(gene_name_list, region_start_list, region_end_list) = associate_island_with_genes.get_gene_feature_region_lists(coords[chrom], opt.region_type, opt.promoter_extension)
			
			tag_start_list = []
			tag_end_list = []
			tag_list = []
			read_file = opt.bedpath + chrom + ".bed";
			f = open(read_file,'r')
			# bed file must be sorted! 
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					tag_list.append(line)
					tag_start_list.append(atoi(sline[1]))
			f.close();
			
			for i in range(0, len(gene_name_list)):
				assert region_start_list[i] <= region_end_list[i]
				o = open(opt.out_file_path + gene_name_list[i] + '.bed', 'w')
				s = bisect.bisect_left(tag_start_list, region_start_list[i])
				e = bisect.bisect_right(tag_start_list, region_end_list[i])
				for j in range(s,e):
					o.write(tag_list[j] + '\n')
				o.close()


if __name__ == "__main__":
	main(sys.argv)

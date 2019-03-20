#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
from GenomeData import *;


plus = re.compile("\+");
minus = re.compile("\-");


def find_upstream_intergenic_region(gene_list, chrom, chrom_length, extension):
	result = []
	for i in range(0, len(gene_list)):
		g = gene_list[i]
		if plus.match(g.strand):
			start = g.txStart - extension
			if i == 0:
				result.append(g.name+'\t'+chrom+'\t-\t0\t'+str(start))
			elif start > (gene_list[i-1].txEnd + extension):
				result.append(g.name+'\t'+chrom+'\t-\t'+str(gene_list[i-1].txEnd + extension)+'\t'+str(start))
		elif minus.match(g.strand):
			start = g.txEnd + extension
			if i == len(gene_list)-1:
				result.append(g.name+'\t'+chrom+'\t+\t'+str(start)+'\t'+str(chrom_length))
			elif start < (gene_list[i+1].txStart-extension):
				result.append(g.name+'\t'+chrom+'\t+\t'+str(start)+'\t'+str(gene_list[i+1].txStart-extension))
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-i", "--input_genes_file", action="store", type="string", dest="genes_file", metavar="<file>", help="input genes file in UCSC known gene format")
	parser.add_option("-e", "--extension", action="store", type="int", dest="extension", metavar="<int>", help="gene region extension in bps")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
		chrom_lengths = species_chrom_lengths[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	coords = UCSC.KnownGenes(opt.genes_file)
	Result = {}
	for chrom in chroms:
		if chrom in coords.keys():
			Result[chrom] = find_upstream_intergenic_region(coords[chrom], chrom, chrom_lengths[chrom], opt.extension)
			print chrom, len(coords[chrom]), len(Result[chrom])
	file = open(opt.out_file, 'w')
	for chrom in Result.keys():
		for g in Result[chrom]:
			file.write(g+'\n')
	file.close()


if __name__ == "__main__":
	main(sys.argv)
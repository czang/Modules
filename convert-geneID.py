#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


"""
This module is used to associate the islands distribution with genes. For each known gene along a genome, it tells: 
1) How many islands are in its promotor region 
2) How many are in its gene body region excluding the promotor region, noted as GeneBody1
3) How many are in its gene body region including the downstream half of the promotor region, noted as GeneBody2

"""


def getNameLists(IDfile):
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if sline[0] != sline[1]:
				refseq[sline[0]] = sline[1]
	return refseq


def convertRefSeqID(genes, RefSeq):
	for gene in genes:
		if RefSeq.has_key(gene['name']):
			gene['name'] = RefSeq[gene['name']]
	return genes


def get_unique_gene_list(coords, RefSeq, species):
	'''Return a unique gene list in one chromosome. Each element in the list consists of gene RefSeq ID, txStart and txEnd. '''
	mm_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX', 'chrY', 'chrM']
	hg_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
	if species == "mm8": 
		chroms = mm_chroms
	elif species == "hg18": 
		chroms = hg_chroms
	else:
		print 'Species Error'
	
	old_name_list = []
	for chrom in chroms: 
		if chrom in coords.keys():
			genes = coords[chrom]
			for g in genes:
				old_name_list.append(g.name)
		else: print 'ERROR!!!';
	old_name_list.sort()
	#print old_name_list
	
	degenerated_old_name_list = []
	degenerated_old_name_list.append(old_name_list[0])
	for i in range(1, len(old_name_list)):
		if old_name_list[i] != old_name_list[i-1]:
			degenerated_old_name_list.append(old_name_list[i])
	#print degenerated_old_name_list
	
	unique_genes = {}
	new_name_list = []
	for chrom in chroms:
		if chrom in coords.keys():
			genes = coords[chrom]
			unique_gene_list = []
			for g in genes:
				#print g.name
				if re.match(g.name, 'BC013184'):
					print g.name; 
				if not RefSeq.has_key(g.name):
					if re.match(g.name, 'BC013184'):
						print '2', g.name
					if not g.name in new_name_list:
						if re.match(g.name, 'BC013184'): print '3', g.name
						new_name_list.append(g.name)
						gene_info = {}
						gene_info['name'] = g.name
						gene_info['strand'] = g.strand
						gene_info['txStart'] = g.txStart
						gene_info['txEnd'] = g.txEnd
						gene_info['cdsStart'] = g.cdsStart
						gene_info['cdsEnd'] = g.cdsEnd
						gene_info['exonCount'] = g.exonCount
						gene_info['exonStarts'] = g.exonStarts
						gene_info['exonEnds'] = g.exonEnds
						unique_gene_list.append(gene_info)
				elif (not RefSeq[g.name] in degenerated_old_name_list) and (not RefSeq[g.name] in new_name_list):
					if re.match(g.name, 'BC013184'):
						print '4', g.name, RefSeq[g.name], new_name_list
					new_name_list.append(RefSeq[g.name])
					gene_info = {}
					gene_info['name'] = RefSeq[g.name]
					gene_info['strand'] = g.strand
					gene_info['txStart'] = g.txStart
					gene_info['txEnd'] = g.txEnd
					gene_info['cdsStart'] = g.cdsStart
					gene_info['cdsEnd'] = g.cdsEnd
					gene_info['exonCount'] = g.exonCount
					gene_info['exonStarts'] = g.exonStarts
					gene_info['exonEnds'] = g.exonEnds
					unique_gene_list.append(gene_info)
			unique_genes[chrom] = unique_gene_list
	return unique_genes


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-k", "--known_genes_file", action="store", type="string", dest="known_genes", metavar="<file>", help="file with known genes")
	parser.add_option("-n", "--gene_names_file", action="store", type="string", dest="gene_names", metavar="<file>", help="file with known gene names converting to ref seq")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	

	coords = UCSC.KnownGenes(opt.known_genes)
	refseq = getNameLists(opt.gene_names)
	print refseq;
	Result = get_unique_gene_list(coords, refseq, opt.species)
	
	f = open(opt.outfile,'w')
	
	for chrom in Result.keys(): 
		for g in Result[chrom]: 
			f.write(g['name'] + '\t' + chrom + '\t' + str(g['strand']) + '\t' + str(g['txStart']) + '\t' + str(g['txEnd']) + '\t' + str(g['cdsStart']) + '\t' + str(g['cdsEnd']) + '\t' + str(g['exonCount']) + '\t' + str(g['exonStarts']) + '\t' + str(g['exonEnds'])+'\n')
	
	f.close()


if __name__ == "__main__":
	main(sys.argv)
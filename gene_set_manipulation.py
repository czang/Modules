#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;



"""
This module is used to compare two sets of gene names and find the same and different genes between them. 
The input are two files, with the RefSeq IDs of the genes in the first column of the first file and in the second column of the second file. 
The output are three files of gene RefSeq lists, which are the same genes and the different genes for each one. 

07/10/09, added the renormalize_column function

"""


def get_gene_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	file = open(gene_file,'r')
	gene_list = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if (len(sline)>c):
				gene_list.append(sline[c])
	file.close()
	return gene_list


def get_float_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[c]))
	file.close()
	return List


def get_int_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(int(atof(sline[c])))
	file.close()
	return List

def renormalize_column(infile, c, const, outfile):
	"""
	c is the 0-based column number 
	"""
	file = open(infile,'r')
	ofile = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			sline[c] = str(atof(sline[c]) * const) 
			line =  '\t '.join(sline) + '\n';
			ofile.write(line);
	file.close()
	ofile.close()
			
			
def get_conversion_table(gene_file, origin, end):
	"""
	origin and end are the 0-based column numbers 
	Return a dictionary
	
	"""
	conversion_table = {};
	file = open(gene_file,'r')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > max(origin, end):
				conversion_table[sline[origin]] = sline[end];
	file.close()
	return conversion_table;


def get_gene_float_dic(expressionfile, colum):
	'''returns a dictionary with geneIDs as keys, expression value as values. colum is the number of colum -1 where the expression data are in the file'''
	file = open(expressionfile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			dic[sline[0]] = atof(sline[colum])
	file.close()
	return dic


def extract_two_columns(gene_file, c1, c2, outfile):
	"""
	c is the 0-based column number 
	"""
	maxi = max (c1, c2);
	mini = min (c1, c2);
	file = open(gene_file,'r')
	ofile = open(outfile, 'w')
	for line in file:
			line = line.strip()
			sline = line.split()
			if len(sline)> maxi:
				outline = sline[c1] + ' ' + sline[c2] + '\n';
			elif len(sline)>mini:
				outline = sline[mini]+ '\n';
			ofile.write(outline);
	ofile.close();
	file.close();
	
	
	
def find_redundant_genes(gene_list):
	"""
	Return the list of names
	"""
	redundance_list=[];
	gene_list.sort();
	if len(gene_list)>=2:
		if gene_list[0] == gene_list[1]:
			redundance_list.append(gene_list[1]);
		for i in range(2, len(gene_list)):
			if gene_list[i] == gene_list[i-1] and gene_list[i-1] != gene_list[i-2]:
				redundance_list.append(gene_list[i]);
	return redundance_list;


def find_unique_genes (gene_list):
	unique_list=[];
	if gene_list !=[]:
		gene_list.sort();
		unique_list.append(gene_list[0]);
		for i in range(1, len(gene_list)):
			if gene_list[i] != gene_list[i-1]:
				unique_list.append(gene_list[i]);
	return unique_list;


def find_unique_only_list (gene_list):
	unique_list=[];
	if gene_list !=[]:
		gene_list.sort();
		if gene_list[1] != gene_list[0]:
			unique_list.append(gene_list[0]);
		for i in range(1, len(gene_list)-1):
			if gene_list[i] != gene_list[i-1] and gene_list[i] != gene_list[i+1]:
				unique_list.append(gene_list[i]);
		if gene_list[len(gene_list)-1] != gene_list[len(gene_list)-2]:
			unique_list.append(gene_list[len(gene_list)-1]);
	return unique_list;


def gene_comparison(List1, List2):
	same = 'shared';
	diff1 = 'only in 1';
	diff2 = 'only in 2';
	
	result_lists = {}
	result_lists[same] = []
	result_lists[diff1] = []
	result_lists[diff2] = []
	for gene in List1:
		if gene in List2:
			result_lists[same].append(gene)
		else:
			result_lists[diff1].append(gene)
	for gene in List2:
		if not gene in result_lists[same]:
			result_lists[diff2].append(gene)
	return result_lists



def get_gene_dictionary(IDfile, c, c1 = 0):
	'''
	'''
	file = open(IDfile,'r')
	refseq = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) > max(c,c1):
				refseq[sline[c1]] = sline[c];
	file.close()
	return refseq


def convertRefSeqID(genes, conversion_table):
	"""
	RefSeq is the conversion table.
	"""
	for gene in genes:
		if conversion_table.has_key(gene['name']):
			gene['name'] = conversion_table[gene['name']]
	return genes


def convertID(list, conversion_table):
	"""
	RefSeq is the conversion table.
	"""
	for gene in list:
		if conversion_table.has_key(gene):
			gene = conversion_table[gene]
	return list


def convert_geneID_in_file(gene_file, geneID_column, conversion_table,  outfile):
	"""
	"""
	file = open(gene_file,'r')
	ofile = open(outfile, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_name = sline[geneID_column];
			#print gene_name;
			if conversion_table.has_key(gene_name):
				sline[geneID_column] = conversion_table[gene_name];	
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: ofile.write(line);
	ofile.close();
	file.close();

def get_UCSCsubset(known_gene_file, specific_gene_list):
	"""
	Return a UCSC Known gene object 
	Get the UCSC file for this  specific_gene_list
	"""
	my_gene_coords ={};
	gene_coords = UCSC.KnownGenes(known_gene_file);
	for chrom in gene_coords.keys():
		my_gene_coords[chrom]=[];
	for my_gene in specific_gene_list:
		for chrom in gene_coords.keys():
			for gene in gene_coords[chrom]:
				if re.match(my_gene, gene.name):
					my_gene_coords[chrom].append(gene);
	return my_gene_coords

def output_UCSCsubset_in_file (known_gene_file, specific_gene_list, outfile):
	"""
	gene_file: is a file of tabular format, must be a UCSC file with name in the 0th column
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be extracted from gene_file and written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	"""
	output_subset_in_file (known_gene_file, 0, specific_gene_list, outfile);
			
	
def output_subset_in_file (gene_file, gene_name_column, specific_gene_list, outfile):
	"""
	gene_file: is a file of tabular format, could be a UCSC file or an expression file
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be extracted from gene_file and written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	This will not take care of redundancy in gene_file
	"""
	file = open(gene_file,'r');
	ofile = open(outfile,'w');
	#print "There are ", len(specific_gene_list), " in your list";
	unique_gene_list = find_unique_genes(specific_gene_list);
	#print len(unique_gene_list)," out of ", len(specific_gene_list), " genes in your list are unique";
	for line in file:
		if not re.match("#", line):
                	line = line.strip();
			sline = line.split();
			if gene_name_column < len(sline):
				for item in unique_gene_list:
					#if re.match(sline[gene_name_column], item):
					if sline[gene_name_column] == item:
						outline = '\t'.join(sline) +'\n';
						ofile.write(outline);
						break; #Once found a match,no need to look further
	file.close();
	ofile.close();


def output_noinsubset_in_file (gene_file, gene_name_column, specific_gene_list, outfile):
	"""
	gene_file: is a file of tabular format, could be a UCSC file or an expression file
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be EXCLUDED from gene_file and NOT written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	This will not take care of redundancy in gene_file
	"""
	file = open(gene_file,'r');
	ofile = open(outfile,'w');
	#print "There are ", len(specific_gene_list), " in your list";
	unique_gene_list = find_unique_genes(specific_gene_list);
	#print len(unique_gene_list)," out of ", len(specific_gene_list), " genes in your list are unique";
	for line in file:
		if not re.match("#", line):
                	line = line.strip();
			sline = line.split();
			if gene_name_column < len(sline):
				if not sline[gene_name_column] in unique_gene_list:
					outline = '\t'.join(sline) +'\n';
					ofile.write(outline);
	file.close();
	ofile.close();


def output_sorted_subset_in_file (gene_file, gene_name_column, specific_gene_list, outfile):
	"""
	07/23/2012: assume both lists are sorted, this one is faster dealing with large data set. 
	gene_file: is a file of tabular format, could be a UCSC file or an expression file
	gene_name_column: specifies the column of the file that is checked for inclusion/exclusion, 0 based
	specific_gene_list; the genes that needs to be extracted from gene_file and written into outfile.
		Even if there are redundancy in specific_gene_list, it should be ok.
	This will not take care of redundancy in gene_file
	"""
	file = open(gene_file,'r');
	ofile = open(outfile,'w');
	#print "There are ", len(specific_gene_list), " in your list";
	#unique_gene_list = find_unique_genes(specific_gene_list);
	#print len(unique_gene_list)," out of ", len(specific_gene_list), " genes in your list are unique";
	i = 0
	L = len(specific_gene_list)
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			if gene_name_column < len(sline):
				target = sline[gene_name_column]
				if i == L:
					break
				elif target == specific_gene_list[i]:
					outline = '\t'.join(sline) +'\n';
					ofile.write(outline);
					i += 1
	file.close();
	ofile.close();

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--genefile1", action="store", type="string", dest="genefile1", metavar="<file>", help="file with gene expression info and RefSeq should be in the second columne")
	parser.add_option("-b", "--genefile2", action="store", type="string", dest="genefile2", metavar="<file>", help="file with island analysis on gene and RefSeq should be in the first columne")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	result = {}
	result = gene_comparison(get_gene_list(opt.genefile1, 1), get_gene_list(opt.genefile2, 0))
	for key in result:
		print key
		for gene in result[key]:
			print gene
		print len(result[key])
		print '\n'



if __name__ == "__main__":
	main(sys.argv)

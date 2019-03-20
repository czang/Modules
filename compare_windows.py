#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

'''
This module is to find the overlapped summary graph windows. Given two summary graph files of the same species, output the number of their overlapped windows.
'''


def compare_lists(List1, List2):
	'''compare two pure number lists with unique elements, return the coexisting number of elements'''
	result = 0
	index1 = 0
	index2 = 0
	List1.sort()
	List2.sort()
	while(index1 < len(List1) and index2 < len(List2)):
		if List1[index1] < List2[index2]:
			index1 += 1
		elif List1[index1] > List2[index2]:
			index2 += 1
		elif List1[index1] == List2[index2]:
			result += 1
			index1 += 1
			index2 += 1
	return result


def compare_bed_graph_lists_by_start(List1, List2):
	'''compare two bed graph lists with unique sorted windows by start positions, return the coexisting number of windows'''
	result = 0
	index1 = 0
	index2 = 0
	while(index1 < len(List1) and index2 < len(List2)):
		if List1[index1].start < List2[index2].start:
			index1 += 1
		elif List1[index1].start > List2[index2].start:
			index2 += 1
		elif List1[index1].start == List2[index2].start:
			result += 1
			index1 += 1
			index2 += 1
	return result


def output_overlapped_bed_by_start(List1, List2):
	Result = []
	outlist1 = [0]*len(List1)
	outlist2 = [0]*len(List2)
	index1 = 0
	index2 = 0
	while(index1 < len(List1) and index2 < len(List2)):
		if List1[index1].start < List2[index2].start:
			index1 += 1
		elif List1[index1].start > List2[index2].start:
			index2 += 1
		elif List1[index1].start == List2[index2].start:
			outlist1[index1] = 1
			outlist2[index2] = 1
			index1 += 1
			index2 += 1
	Result.append(outlist1)
	Result.append(outlist2)
	return Result


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summarygraphfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="summary graph file 1")
	parser.add_option("-b", "--summarygraphfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="summary graph file 2")
	parser.add_option("-t", "--output_type", action="store", type="int", dest="outtype", help="0 for a number, 1 for two overlapped files", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	summary1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0);
	summary2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0);
	
	if opt.outtype == 0:
		Result = 0
		for chrom in chroms:
			if chrom in summary1.keys() and chrom in summary2.keys():
				Result += compare_bed_graph_lists_by_start(summary1[chrom], summary2[chrom])
		print Result
	
	elif opt.outtype ==1:
		f = open(opt.bedfile1+'_'+opt.bedfile2+'_overlapped','w')
		g = open(opt.bedfile2+'_'+opt.bedfile1+'_overlapped','w')
		for chrom in chroms:
			if chrom in summary1.keys() and chrom in summary2.keys():
				ref = output_overlapped_bed_by_start(summary1[chrom], summary2[chrom])
				for i in range(0, len(summary1[chrom])):
					if ref[0][i] == 1:
						f.write(chrom + '\t' + str(summary1[chrom][i].start) + '\t' + str(summary1[chrom][i].end) + '\t' + str(summary1[chrom][i].value) + '\n')
				for i in range(0, len(summary2[chrom])):
					if ref[1][i] == 1:
						g.write(chrom + '\t' + str(summary2[chrom][i].start) + '\t' + str(summary2[chrom][i].end) + '\t' + str(summary2[chrom][i].value) + '\n')
		f.close()
		g.close()
	
	else:
		print "Output type error";
		sys.exit(1);


if __name__ == "__main__":
	main(sys.argv)

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

## get BED module
import BED
from GenomeData import *
import separate_by_chrom

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");


def make_summary_graph_dic(chrom, file):
	dic = {}
	f = open(file,'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if sline[0] == chrom:
				dic[atoi(sline[1])] = atof(sline[3])
	return dic


def window_tag_density(dic, chrom_length, window_size):
	total_tag = 0.0
	for position in dic.keys():
		total_tag += dic[position]
	return total_tag / float(chrom_length) * window_size


def correlation2(dic1, dic2, chrom_length, dx, r):
	average1 = window_tag_density(dic1, chrom_length, dx)
	average2 = window_tag_density(dic2, chrom_length, dx)
	keylist1 = dic1.keys()
	keylist1.sort()
	keylist2 = dic2.keys()
	keylist2.sort()
	x = 0
	i = 0
	j = min(len(keylist2)-1, bisect.bisect_left(keylist2,int(r)))
	total = 0.0
	while x < (chrom_length-r):
		if x == keylist1[i]:
			Tx = dic1[x]
			if i < (len(keylist1)-1):
				i += 1
		elif x < keylist1[i]:
			Tx = 0
		elif x > keylist1[i]:
			Tx = 0
			if i < (len(keylist1)-1):
				i += 1
		y = int(x + r)
		if y == keylist2[j]:
			Tr = dic2[y]
			if j < (len(keylist2)-1):
				j += 1
		elif y < keylist2[j]:
			Tr = 0
		elif y > keylist2[j]:
			Tr = 0
			if j < len(keylist2)-1:
				j += 1
		total += (float(Tx))*(float(Tr))
		x += dx
	return total/chrom_length*dx


def correlation_function2(dic1, dic2, chrom_length, step, dx):
	result = {}
	r = 0
	while r < min(chrom_length, step*11):
		result[r] = correlation2(dic1, dic2, chrom_length, dx, r)
		r += step
	return result


def generate_all_functions_2(extension1, extension2, chroms, dx, dr, file_name, chrom_lengths):
	for chrom in chroms:
		dic1 = make_summary_graph_dic(chrom, chrom+extension1)
		dic2 = make_summary_graph_dic(chrom, chrom+extension2)
		chrom_length = chrom_lengths[chrom]
		if len(dic1.keys()) > 0 and len(dic2.keys()) > 0:
			result_dic = correlation_function2(dic1, dic2, chrom_length, dr, dx)
			keylist = result_dic.keys()
			keylist.sort()
			f = open(chrom+'_'+file_name, 'w')
			for i in keylist:
				f.write(str(i)+'\t'+str(result_dic[i])+'\n')
			f.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summary_graph_file1", action="store", type="string", dest="bedfile1", metavar="<file>", help="summary graph file 1 in bed format")
	parser.add_option("-b", "--summary_graph_file2", action="store", type="string", dest="bedfile2", metavar="<file>", help="summary graph file 2 in bed format")
	parser.add_option("-i", "--windows_size", action="store", type="int", dest="window_size", metavar="<int>", help="window size in summary graph file")
	parser.add_option("-d", "--data_resolution", action="store", type="int", dest="step", metavar="<int>", help="distance between data points, must be integer times of window size")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file extension")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
		chrom_lengths = species_chrom_lengths[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separate_by_chrom.separateByChrom(chroms, opt.bedfile1, '.bed1')
	separate_by_chrom.separateByChrom(chroms, opt.bedfile2, '.bed2')
	generate_all_functions_2('.bed1', '.bed2', ['chr1'], opt.window_size, opt.step, opt.out_file, chrom_lengths)
	separate_by_chrom.cleanup(chroms, '.bed1')
	separate_by_chrom.cleanup(chroms, '.bed2')


if __name__ == "__main__":
	main(sys.argv)

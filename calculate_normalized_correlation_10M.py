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


def correlation1(dic, chrom_length, dx, r):
	average = window_tag_density(dic, chrom_length, dx)
	keylist = dic.keys()
	keylist.sort()
	x = 0
	i = 0
	j = min(len(keylist)-1, bisect.bisect_left(keylist,int(r)))
	total = 0.0
	k = 0
	while x < (chrom_length-r):
		if x == keylist[i]:
			Tx = dic[x]
			if i < (len(keylist)-1):
				i += 1
		elif x < keylist[i]:
			Tx = 0
		elif x > keylist[i]:
			Tx = 0
			if i < (len(keylist)-1):
				i += 1
		xr = int(x + r)
		if xr == keylist[j]:
			Tr = dic[xr]
			if j < (len(keylist)-1):
				j += 1
		elif xr < keylist[j]:
			Tr = 0
		elif xr > keylist[j]:
			Tr = 0
			if j < len(keylist)-1:
				j += 1
		total += (float(Tx) - average)*(float(Tr) - average)
		k += 1
		x += dx
	return total/k


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


def correlation_function(dic, chrom_length, step, dx):
	result = {}
	normalization = correlation1(dic, chrom_length, dx, 0)
	print normalization
	r = 100000
	while r <= min(chrom_length, 10000000):
		result[r] = correlation1(dic, chrom_length, dx, r) / normalization
		r += step
	return result


def generate_all_functions(chroms, dx, dr, file_name, chrom_lengths):
	for chrom in chroms:
		dic = make_summary_graph_dic(chrom, chrom+".bed")
		chrom_length = chrom_lengths[chrom]
		if len(dic.keys()) > 0:
			print chrom
			result_dic = correlation_function(dic, chrom_length, dr, dx)
			keylist = result_dic.keys()
			keylist.sort()
			f = open(chrom+file_name, 'w')
			for i in keylist:
				f.write(str(i)+'\t'+str(result_dic[i])+'\n')
			f.close()


def cleanup(chroms):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        try:
            if os.remove('%s' %
                         (bed_file)): raise
        except: sys.stderr.write("clean up failed\n");


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--summary_graph_file", action="store", type="string", dest="bedfile", metavar="<file>", help="summary graph file in bed format")
	parser.add_option("-i", "--windows_size", action="store", type="int", dest="window_size", metavar="<int>", help="window size in summary graph file")
	parser.add_option("-d", "--data_resolution", action="store", type="int", dest="step", metavar="<int>", help="distance between data points, must be integer times of window size")
	#parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="filtered summary graph file")
	
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
	
	separate_by_chrom.separateByChrom(chroms, opt.bedfile, ".bed")
	generate_all_functions(chroms, opt.window_size, opt.step, "_correlation.txt", chrom_lengths)
	#filter_tags_by_islands(chroms, islands, opt.fragment_size)
	#final_output_file = opt.out_file;
	#final_output_file = combineAllGraphFiles(chroms, final_output_file);
	cleanup(chroms);


if __name__ == "__main__":
	main(sys.argv)

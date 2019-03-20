#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from GenomeData import *
import Utility

def is_list_sorted(List):
        """
        Check if sorted in ascending order.
        input is a list of values.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(List)-1):
                if List[index] > List[index+1]:
                        sorted = 0;
        return sorted;

def detect_dip_from_list(List, t1 = 0.1):
	d1 = derivative(List)
	d2 = derivative(d1)
	index_result = []
	for i in range(2, len(List)-2):
		i1 = i - 1
		i2 = i1 - 1
		m = float(List[i] + List[i1] + List[i2]) / 3
		if d1[i1] <= t1 * m and d2[i2] > 0:
			index_result.append(i)
	return index_result

def detect_peak_from_list(List, t1 = 0.1):
	d1 = derivative(List)
	d2 = derivative(d1)
	index_result = []
	for i in range(2, len(List)-2):
		i1 = i - 1
		i2 = i1 - 1
		m = float(List[i] + List[i1] + List[i2]) / 3
		if d1[i1] <= t1 * m and d2[i2] < 0:
			index_result.append(i)
	return index_result

def get_sublist_index(List, start, end):
	if is_list_sorted(List) == 0:
		List.sort()
	i = bisect_left(List, start)
	j = bisect_right(List, end)
	return i, j

def wiglists2datalist(List1,List2, step=10):
	assert len(List1) == len(List2)
	result = []
	result.append(List2[0])
	flag = List1[0]
	i = 1
	while i < len(List1):
		if List1[i] == flag + step: 
			result.append(List2[i])
			flag = List1[i]
			i += 1
		else:
			while List1[i] >= flag + step:
				result.append(0)
				flag += step
			i += 1
	assert len(result) * (step-1) == List1[-1] - List[1]
	return result

def derivative(List):
	result = []
	for i in range(1,len(List)-1):
		result.append(float(List[i+1] - List[i-1])/2)
	assert len(result) + 2 == len(List)
	return result

def peak_5(List):
	assert len(List) == 5
	if List[2] > List[1] and List[2] > List[3]:
		if List[1] > List[0] and List[3] > List[4]:
			return 1
		else:
			return -1
	else:
		return -1

def dip_5(List):
	assert len(List) == 5
	if List[2] < List[1] and List[2] < List[3]:
		if List[1] < List[0] and List[3] < List[4]:
			return 1
		else:
			return -1
	else:
		return -1


def wig_to_lists(file):
	positionlist = []
	countlist = []
	f = open(file,'r')
	for line in f:
		if not (re.match("track", line) or re.match("variableStep",line)):
			line = line.strip()
			sline = line.split()
			positionlist.append(atoi(sline[0]))
			countlist.append(atoi(sline[1]))
	f.close()
	return positionlist, countlist


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-g", "--wiggle", action="store", type="string", dest="wig", metavar="<file>", help="wig files")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	o = open(opt.out_file, 'w')
	for chrom in chroms:
		if Utility.fileExists(chrom + opt.wig):
			print 'reading', chrom, 'wig ...'
			(positionlist, countlist) = wig_to_lists(chrom + opt.wig)
			for i in range(0,len(countlist)-5):
				if positionlist[i+4] - positionlist[i] < 50: 
					List = countlist[i:i+5]
					if peak_5(List) == 1: 
						start = positionlist[i+2] - 1
						end = start + 1
						score = List[2]
						o.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)

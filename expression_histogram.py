#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *

Dir = os.getcwd();


def get_float_data_list(gene_file, c):
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
			gene_list.append(atof(sline[c]))
	file.close()
	return gene_list


def histogram(List, step):
	List.sort()
	histogramlist = {}
	i = 0
	down = 0
	up = down + step
	while i < (len(List)):
		k = 0
		while i < (len(List)) and List[i] >= down and List[i] < up:
			k+=1
			i+=1
		histogramlist[up] = float(k)/float(len(List))*100
		down += step
		up += step
	return histogramlist


def histogram_change_total(List, step, total):
	List.sort()
	histogramlist = {}
	i = 0
	down = 0
	up = down + step
	while i < (len(List)):
		k = 0
		while i < (len(List)) and List[i] >= down and List[i] < up:
			k+=1
			i+=1
		histogramlist[up] = float(k)/float(total)*100
		down += step
		up += step
	return histogramlist


def cumulative_histogram(List, step):
	List.sort()
	result = []
	k = 0
	while bisect.bisect_left(List, k) < len(List):
		result.append(float(len(List)-bisect.bisect_left(List, k))/float(len(List))*100)
		k +=step
	return result


def probability_distribution(List, bin):
	List.sort()
	result = []
	down = 0.0
	up = float(bin)
	leftup = bisect.bisect_left(List, up)
	while leftup <= len(List):
		number = leftup - bisect.bisect_left(List, down)
		#result.append(number)
		result.append(float(number)/len(List)*100)
		down += bin
		up += bin
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-c", "--colum", action="store", type="int", dest="colum", metavar="<int>", help="colum number for data, starts from 0")
	parser.add_option("-w", "--bin", action="store", type="float", dest="step", metavar="<float>", help="bin size on x axis of the histogram")
	parser.add_option("-t", "--type", action="store", type="int", dest="resulttype", metavar="<int>", help="1 for integral curve, 0 for histogram")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	datalist = get_float_data_list(opt.infile, opt.colum)
	#genelist = get_list(opt.infile, 0)
	if opt.resulttype == 1:
		plotdata = cumulative_histogram(datalist, opt.step)
		f = open(opt.outfile,'w')
		for i in range(0, len(plotdata)):
			f.write(str(i*opt.step)+'\t'+str(plotdata[i])+'\n')
		f.close()
	elif opt.resulttype == 0:
		Histogram = histogram(datalist, opt.step)
		f = open(opt.outfile,'w')
		for i in Histogram.keys():
			if Histogram[i] != 0:
				f.write(str(i)+'\t'+str(Histogram[i])+'\n')
		f.close()
	elif opt.resulttype == -1:
		Histogram = probability_distribution(datalist, opt.step)
		f = open(opt.outfile,'w')
		for i in range(0, len(Histogram)):
			f.write(str(i*opt.step)+'\t'+str(Histogram[i])+'\n')
		f.close()
	else:
		print "input error!"


if __name__ == "__main__":
	main(sys.argv)
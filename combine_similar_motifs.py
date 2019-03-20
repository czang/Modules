#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-s", "--scoreindex", action="store", type="int", dest="score_index", metavar="<int>", help="column for score, starts from 0")
	parser.add_option("-c", "--clusterindex", action="store", type="int", dest="cluster_index", metavar="<int>", help="column for cluster, starts from 0")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	Total = 245
	score_index = opt.score_index
	cluster_index = opt.cluster_index
	List = [''] * Total
	score_list = [-0.1] * Total
	inputfile = open(opt.infile,'r')
	for line in inputfile:
		if not (re.match("#", line) or re.match("track", line) or re.match("variableStep", line)):
			line = line.strip()
			sline = line.split()
			assert len(sline) >= max(score_index, cluster_index)
			i = atoi(sline[cluster_index]) - 1
			assert i < Total
			score = atof(sline[score_index])
			if score > score_list[i]:
				List[i] = sline[0]
				score_list[i] = score
	inputfile.close()
	outfile = open(opt.out_file, 'w')
	for i in range(0, Total):
		outfile.write(List[i] + '\t' + str(score_list[i]) + '\n')

	
if __name__ == "__main__":
	main(sys.argv)

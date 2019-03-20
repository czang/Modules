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
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""
ch = {}
ch['1']='chr1'
ch['2']='chr2'
ch['3']='chr3'
ch['4']='chr4'
ch['5']='chr5'
ch['6']='chr6'
ch['7']='chr7'
ch['8']='chr8'
ch['9']='chr9'
ch['10']='chr10'
ch['11']='chr11'
ch['12']='chr12'
ch['13']='chr13'
ch['14']='chr14'
ch['15']='chr15'
ch['16']='chr16'
ch['17']='chr17'
ch['18']='chr18'
ch['19']='chr19'
ch['20']='chr20'
ch['21']='chr21'
ch['22']='chr22'
ch['X']='chrX'
ch['Y']='chrY'
ch['M']='chrM'


def chromname(number):
	if number in ch.keys():
		return ch[number]
	else:
		return number


def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number


def get_formated_genes(input_file, name, chrom, strand, start, end, output_file):
	"""
	c is the 0-based column number 
	Return a list of names
	
	"""
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	#gene_list = []
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= max(name, chrom, strand, start, end):
				if strand > 0:
					outfile.write(sline[name - 1] + '\t' + chromname(sline[chrom - 1]) + '\t' + strandsign(sline[strand - 1]) + '\t' + sline[start - 1] + '\t' + sline[end - 1] + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\n')
				else:
					outfile.write(sline[name - 1] + '\t' + chromname(sline[chrom - 1]) + '\t+\t' + sline[start - 1] + '\t' + sline[end - 1] + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\t' + '0' + '\n')
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-a", "--genename", action="store", type="int", dest="name", metavar="<int>", help="colum number for gene names")
	parser.add_option("-b", "--chromosome", action="store", type="int", dest="chrom", metavar="<int>", help="colum number for chromosome")
	parser.add_option("-c", "--strand", action="store", type="int", dest="strand", metavar="<int>", help="colum number for gene strand, all + if 0")
	parser.add_option("-d", "--txStart", action="store", type="int", dest="start", metavar="<int>", help="colum number for transcription start position")
	parser.add_option("-e", "--txEnd", action="store", type="int", dest="end", metavar="<int>", help="colum number for transcription end position")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
	
	get_formated_genes(opt.input_file, opt.name, opt.chrom, opt.strand, opt.start, opt.end, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)

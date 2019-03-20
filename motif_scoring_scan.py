#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import SeparateByChrom
import Utility
from GenomeData import *


genomefasta_dir = '/cluster2/czang/genomes/humanhg18/rawgenome/'

DNase_BG_0 = [0.28031827806535914, 0.2193352207867304, 0.21958461937786461, 0.28076188177004585]

def pwmfile2list(file):
	f = open(file, 'r')
	m = []
	for line in f:
		if not re.search("A",line):
			line = line.strip()
			sline = line.split()
			assert len(sline) == 4
			List = []
			for item in sline:
				List.append(atof(item))
			m.append(List)
	f.close()
	return m


def nt2number(read):
	temp0 = read.replace('A','0')
	temp1 = temp0.replace('C','1')
	temp2 = temp1.replace('G','2')
	temp3 = temp2.replace('T','3')
	temp0 = temp3.replace('a','0')
	temp1 = temp0.replace('c','1')
	temp2 = temp1.replace('g','2')
	temp3 = temp2.replace('t','3')
	return temp3


def nt2number_reversecomplement(read):
	temp0 = read.replace('A','3')
	temp1 = temp0.replace('C','2')
	temp2 = temp1.replace('G','1')
	temp3 = temp2.replace('T','0')
	temp0 = temp3.replace('a','3')
	temp1 = temp0.replace('c','2')
	temp2 = temp1.replace('g','1')
	temp3 = temp2.replace('t','0')
	return temp3[::-1]


def pwmlist2log(List):
	m = []
	for item in List:
		p = []
		for n in item:
			p.append(log(n))
		m.append(p)
	return m


def list2log(List):
	bg = []
	for item in List:
		bg.append(log(item))
	return bg


def allacgt(String):
	nt = 'acgtACGT'
	flag = 1
	i = 0
	while i < len(String):
		if not String[i] in nt:
			flag = 0
			break
		else:
			i += 1
	return flag


def all0123(String):
	nt = '0123'
	flag = 1
	i = 0
	while i < len(String):
		if not String[i] in nt:
			flag = 0
			break
		else:
			i += 1
	return flag


'''from a 0123 sequence with length = motif, calculate log likelihood ratio'''
def logLikelihoodRatio(seq, logpwmlist, logBG0list):
	'''seq can only contain 0,1,2,3'''
	length = len(seq)
	assert length == len(logpwmlist)
	llr = 0.0
	for i in range(0, length):
		'''i is index, n is label 0,1,2,3'''
		n = atoi(seq[i])
		ratio = logpwmlist[i][n] - logBG0list[n]
		llr += ratio
	return llr


'''scan one sequence and get the log2 sum likelihood ratio score'''
def logSumLikelihoodRatio(sequence, logpwmlist, logBG0list):
	size = len(logpwmlist)
	length = len(sequence)
	if length >= size: 
		numsequence = nt2number(sequence)
		reverse_comp_num = nt2number_reversecomplement(sequence)
		sum_score = 0.0
		for i in range(0, length - size + 1):
			seq = numsequence[i:i+size]
			'''check whether the sequence contains only ACGT(acgt)'''
			if all0123(seq) == 1:
				llr = logLikelihoodRatio(seq, logpwmlist, logBG0list)
				lr = exp(llr)
				sum_score += lr
			seq = reverse_comp_num[i:i+size]
			'''check whether the sequence contains only ACGT(acgt)'''
			if all0123(seq) == 1:
				llr = logLikelihoodRatio(seq, logpwmlist, logBG0list)
				lr = exp(llr)
				sum_score += lr
		return log(sum_score, 2)
	else:
		return 0.0


def calculate_mm0_bg(bedfile, chroms):
	SeparateByChrom.separateByChrom(chroms, bedfile, '.bed')
	for chrom in chroms:
		chromfasta = genomefasta_dir + chrom + '.fa'
		tmpbed = chrom + '.bed'
		tmpfa = chrom + '.tmp.fa'
		try:
			if os.system('/cluster2/czang/Software/BEDTools/bin/fastaFromBed %s %s %s %s %s %s' % ('-db', chromfasta, '-bed', tmpbed, '-fo', tmpfa)): raise
		except: sys.stderr.write("get fasta sequence failed\n")
	fafile = 'temp.fa'
	SeparateByChrom.combineAllGraphFiles(chroms, '.tmp.fa', fafile)
	SeparateByChrom.cleanup(chroms, '.tmp.fa')
	SeparateByChrom.cleanup(chroms, '.bed')
	
	List = [0]*4
	if Utility.fileExists(fafile):
		f = open(fafile, 'r')
		for line in f:
			if not re.match('>', line): 
				line = line.strip()
				List[0] += line.count('A')
				List[0] += line.count('a')
				List[1] += line.count('C')
				List[1] += line.count('c')
				List[2] += line.count('G')
				List[2] += line.count('g')
				List[3] += line.count('T')
				List[3] += line.count('t')
		s = atof(sum(List))
		for i in range(0, len(List)):
			List[i] = List[i]/s
	return List


def scan_fasta(fastafile, pwmlist, BG0list, outfile):
	logpwmlist = pwmlist2log(pwmlist)
	logBG0list = list2log(BG0list)
	o = open(outfile, 'w')
	f = open(fastafile, 'r')
	for line in f:
		if  re.match('>', line):
			line = line.strip()
			line = line.strip('>')
			line = line.replace(':','\t')
			line = line.replace('-','\t')
			line = line.replace('_','\t')
			o.write(line + '\t')
		else:
			line = line.strip()
			score = logSumLikelihoodRatio(line, logpwmlist, logBG0list)
			o.write(str(score) + '\n')
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, e.g. mm8, hg18", metavar="<str>")
	parser.add_option("-p", "--file1", action="store", type="string", dest="file1", metavar="<file>", help="file 1")
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="regions to be scanned")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)

	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	m = pwmfile2list(opt.file1)
	print m
	
	BGlist = calculate_mm0_bg(opt.bedfile, chroms)
	print BGlist
	
	scan_fasta('temp.fa', m, BGlist, opt.outfile)
	

if __name__ == "__main__":
	main(sys.argv)
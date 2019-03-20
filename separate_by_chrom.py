#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

Dir = os.getcwd();
grep = "/usr/bin/grep";
cat = "/bin/cat";

plus = re.compile("\+");
minus = re.compile("\-");


def separateByChrom(chroms, file, extension):
	for chrom in chroms:
		match = chrom + "[[:space:]]";
		tmpFile = chrom + extension;
		try:
			if os.system('%s %s %s > %s' % (grep, match, file, tmpFile)): raise
		except: sys.stderr.write("grep failed\n");


def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def combineAllGraphFiles(chroms, extension, final_out):
	"""
		Combine the seperately processed chromosomes, return the output file name
	"""
	outfile = open(final_out,'w');
	outfile.close();
    
	for chrom in chroms:
		file = chrom + extension;
		try:
			if os.system('%s %s >> %s' %(cat, file, final_out)): raise
		except: sys.stderr.write("cat failed\n")
	return final_out


def cleanup(chroms, extension):
	for chrom in chroms:
		file = chrom + extension;
		try:
			if os.remove('%s' %(file)): raise
		except: sys.stderr.write("clean up failed\n");


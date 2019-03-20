#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def multiplication(column, in_file, normalization_factor, out_file):
	"""
	column is the 0-based column number 
	"""
	assert ( normalization_factor != 0);
	file = open(in_file,'r');
	ofile = open(out_file, 'w')
	number_of_column=0;
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			sline[column] = str (atof(sline[column])/normalization_factor);
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: 
			ofile.write(line);
	ofile.close();
	file.close();
	

def normalization(list, normalization_factor):
	return [elem/normalization_factor for elem in list];

def addition(column, in_file, addition_factor, out_file):
	"""
	column is the 0-based column number 
	"""
	file = open(in_file,'r');
	ofile = open(out_file, 'w')
	number_of_column=0;
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			sline[column] = str ( atof(sline[column])+addition_factor );
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: 
			ofile.write(line);
	ofile.close();
	file.close();
	
def add_two_columns(column1, column2, in_file, out_file):
	"""
	column is the 0-based column number 
	This one add column2 onto column1
	
	"""
	file = open(in_file,'r');
	ofile = open(out_file, 'w')
	number_of_column=0;
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			sline[column1] = str (atof(sline[column1]) + atof(sline[column2]));
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: 
			ofile.write(line);
	ofile.close();
	file.close();

def multiply_two_columns(column1, column2, in_file, out_file):
	"""
	column is the 0-based column number 
	This one multiplies column2 onto column1
	
	"""
	file = open(in_file,'r');
	ofile = open(out_file, 'w')
	number_of_column=0;
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			sline[column1] = str (atof(sline[column1]) * atof(sline[column2]));
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: 
			ofile.write(line);
	ofile.close();
	file.close();	

def remove(column, in_file, out_file):
	"""
	column is the 0-based column number 
	This one delete the specified column
	
	"""
	file = open(in_file,'r');
	ofile = open(out_file, 'w')
	for line in file:
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			sline.pop(column); 
			outline =  '\t '.join(sline) + '\n';
			ofile.write(outline);
		else: 
			ofile.write(line);
	ofile.close();
	file.close();		


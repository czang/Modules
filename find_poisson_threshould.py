#!/usr/bin/python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

def factorial(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;

def poisson(i, average):
	return exp(-average) * pow(average,i) / factorial(i);

def find_threshold(pvalue, average):
	value = 1;
	index = 0;
	value -= poisson(index, average);
	while value > pvalue:
		index += 1;
		value -= poisson(index, average);	
	return index;


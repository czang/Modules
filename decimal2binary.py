#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

def Denary2Binary(n):
	'''convert denary integer n to binary string bStr'''
	bStr = ''
	if n < 0: raise ValueError, "must be a positive integer"
	if n == 0: return '0'
	while n > 0:
		bStr = str(n % 2) + bStr
		n = n >> 1
	return bStr


def Denary2Binary_38(n):
	Str = Denary2Binary(n)
	if len(Str)<38:
		for i in range(0, 38-len(Str)):
			Str = '0'+Str
	return(Str)


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--instring", action="store", type="string", dest="instring", metavar="<str>", help="input")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	string = opt.instring.split('_')
	print 'mv', opt.instring, string[0]+'_'+string[1]+'_'+Denary2Binary(atoi(string[2]))+'_'+string[3]

if __name__ == "__main__":
	main(sys.argv)
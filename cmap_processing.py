#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Author: Chongzhi Zang
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# czang@jimmy.harvard.edu).

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator


## from a input data list of id, concentration, drug name from one cell type. output 1019 drug samples: id, drug. 


def item_in_list(a, List):
	for item in List:
		if re.match(a, item):
			return True
	return False


def get_list(infile):
	f = open(infile, 'r')
	List = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			List.append(line)
	f.close()
	return List

'''
def process(infile, druglist, outfile):
	f = open(infile,'r')
	id_list = []
	concentration_list = []
	alldrug_list = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			drug = ' '.join(sline[2:])
			#print drug
			#if item_in_list(drug, druglist):
			if drug in druglist:
				id_list.append(sline[0])
				concentration_list.append(atof(sline[1]))
				alldrug_list.append(drug)
	f.close()
	o = open(outfile, 'w')
	#current_drug = alldrug_list[0]
	#current_con = concentration_list[0]
	for i in range(0, len(id_list)):
		o.write(id_list[i]+'\t'+alldrug_list[i]+'\n')
	
	#for batch in datadic.keys():
	#	o.write(datadic[batch]+'\n')
	#f.close()
	o.close()
'''

def process(infile, druglist, outfile):
	f = open(infile,'r')
	id_dic = {}
	con_dic = {}
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			drug = ' '.join(sline[2:])
			#print drug
			#if item_in_list(drug, druglist):
			if drug in druglist:
				if drug in id_dic.keys():
					if atof(sline[1]) > con_dic[drug]:
						id_dic[drug] = sline[0]
						con_dic[drug] = atof(sline[1])
				else:
					id_dic[drug] = sline[0]
					con_dic[drug] = atof(sline[1])
				#id_list.append(sline[0])
				#concentration_list.append(atof(sline[1]))
				#alldrug_list.append(drug)
	f.close()
	o = open(outfile, 'w')
	#current_drug = alldrug_list[0]
	#current_con = concentration_list[0]
	for drug in id_dic.keys():
		o.write(id_dic[drug]+'\t'+drug+'\n')
	
	#for batch in datadic.keys():
	#	o.write(datadic[batch]+'\n')
	#f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-a", "--all_input_file", action="store", type="string",
                      dest="in_file", help="island bed file", metavar="<file>")
	parser.add_option("-b", "--druglist_file", action="store", type="string",
                      dest="drug_file", help="island bed file", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	druglist = get_list(opt.drug_file)
	print len(druglist)
	process(opt.in_file, druglist, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)

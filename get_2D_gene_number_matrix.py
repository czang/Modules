#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def read_data_file(datafile):
	'''returns a list of dictionaries. Each dictionary has three keys: name, x and y.'''
	file = open(datafile,'r')
	result = []
	for line in file:
		item = {}
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			item['name'] = str(sline[0])
			item['x'] = atof(sline[1])
			item['y'] = atof(sline[2])
			result.append(item)
	file.close()
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--datafile", action="store", type="string", dest="inputfile", metavar="<file>", help="input file with complete data, gene, x-value, y-value")
	parser.add_option("-n", "--intervalnumber", action="store", type="int", dest="groupnumber", metavar="<int>", help="number of intevals on one side")
	parser.add_option("-x", "--xmax", action="store", type="float", dest="xmax", metavar="<float>", help="upper limit of x value")
	parser.add_option("-y", "--ymax", action="store", type="float", dest="ymax", metavar="<float>", help="upper limit of y value")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output matrix file name")
	parser.add_option("-d", "--output_data_file", action="store", type="string", dest="out_data", metavar="<file>", help="output data file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	xstep = opt.xmax / opt.groupnumber
	ystep = opt.ymax / opt.groupnumber
	
	ylist = []
	for i in range(0, opt.groupnumber):
		ylist.append(-1 * opt.ymax + i * ystep)
	for i in range(0, opt.groupnumber):
		ylist.append(i * ystep)
	ylist.append(opt.ymax)
	print "y:", ylist
	
	xlist = []
	for i in range(0, opt.groupnumber):
		xlist.append(-1 * opt.xmax + i * xstep)
	for i in range(0, opt.groupnumber):
		xlist.append(i * xstep)
	xlist.append(opt.xmax)
	print "x:", xlist
	
	all_data = read_data_file(opt.inputfile)
	
	#y_sorted = sorted(all_data, key=operator.itemgetter('y'))
	
	y_data = []
	for i in range(0, (len(ylist) + 1)):
		List = []
		y_data.append(List)
	
	y_0_data = []
	
	#c = 0
	for item in all_data:
		if item['y'] < -1E-20:
			index = bisect.bisect_right(ylist, item['y'])
			List = y_data[index]
			List.append(item)
			#for i in range(0,len(y_data)):
				#print len(y_data[i])
		elif item['y'] > 1E-20:
			index = bisect.bisect_left(ylist, item['y'])
			y_data[index].append(item)
			#print index
			#c +=1
		else:
			y_0_data.append(item)
	
	#print c
	#print len(y_data)
	#for i in range(0,len(y_data)):
	#	print len(y_data[i])
	
	x_data = []
	for i in range(0, (len(ylist) + 1)):
		List1 = []
		for j in range(0, (len(xlist) + 1)):
			List2 = []
			List1.append(List2)
		x_data.append(List1)
	
	x_0_data = []
	for i in range(0, (len(ylist) + 1)):
		List = []
		x_0_data.append(List)
	#print len(x_0_data)
	
	for i in range(0, len(y_data)):
		for item in y_data[i]:
			if item['x'] < -1E-20:
				index = bisect.bisect_right(xlist,item['x'])
				x_data[i][index].append(item)
			elif item['x'] > 1E-20:
				index = bisect.bisect_left(xlist,item['x'])
				x_data[i][index].append(item)
			else:
				x_0_data[i].append(item)
		#for j in range(0, len(x_data[i])):
			#print len(x_data[i][j])
	
	# output data
	f = open(opt.out_file, 'w')
	for i in range(0, len(x_data)):
		index = -1 * i -1
		#index = i
		for j in range(0, len(x_data[index])):
			f.write(str(len(x_data[index][j])) + '\t')
		f.write('\n')
	f.close()
	g = open(opt.out_data, 'w')
	for i in range(0, len(x_data)):
		for j in range(0, len(x_data[i])):
			g.write(str(j-opt.xmax) + '\t' + str(i-opt.ymax) + '\t' + str(len(x_data[i][j])) + '\n')
	g.close()
	
	f = open(opt.out_file, 'r')
	for line in f:
		print line
	f.close()



if __name__ == "__main__":
	main(sys.argv)
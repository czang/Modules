#!/usr/bin/env python

def horizontal2vertical(infile, outfile):
	f = open(infile, 'r')
	o = open(outfile, 'w')
	List = []
	for line in f:
		line = line.strip()
		List.append(line.split())
	for i in range(0, len(List[0])):
		o.write(List[0][i]+'\t'+List[1][i]+'\n')
	o.close()
	f.close()
#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator

import BED;
import UCSC;
import bisect;
import GenomeData;

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";
plus = re.compile("\+");
minus = re.compile("\-");

# Many of the modules in this code needs a sorted input! #


def breakUpStrands(bed_list):
	"""
	input: a list of bed6 object 
	output: two lists of bed6 objects, one for plus strand, one for minus strand. 
	"""

	plus_bed_list = [];
	minus_bed_list = [];
	for b in bed_list:
		if plus.match(b.strand):
			plus_bed_list.append(b);
		elif minus.match(b.strand):
			minus_bed_list.append(b);
	return (plus_bed_list, minus_bed_list)



def find_read_copy_distribution(sorted_bed_list):
	"""
	Input:  a list of sorted bed6 objects. Already assumed that 
	the tags are from one chromosome and in one direction. 
	
	Output: the histogram of the tag copies
	
	"""
	unique_tag_histogram = [0] * 100;
	if (len(sorted_bed_list) != 0):
		total_number_tags = len(sorted_bed_list);
		#sorted_bed_list.sort(key=operator.attrgetter('start'));
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (len(unique_tag_histogram)-1)<current_count:
					unique_tag_histogram +=[0]*(current_count-(len(unique_tag_histogram)-1));
				unique_tag_histogram[current_count] += 1;
				current_value = item.start;
				current_count = 1; #reset
			else:
				current_count += 1;
		#last read	
		if (len(unique_tag_histogram)-1)<current_count:
			unique_tag_histogram +=[0]*(current_count-(len(unique_tag_histogram)-1));
		unique_tag_histogram[current_count] += 1; 
	return unique_tag_histogram;
		
def find_multi_copy_reads(sorted_bed_list, threshold):
	"""
	Input:  a list of sorted bed6 objects. Already assumed that 
	the tags are from one chromosome and in one direction. 
		the threshold for read copy
	Output: the list of BED6 reads with copy number above or equal threshold.
	"""
	multiple_copy_read_list=[];
	
	if (len(sorted_bed_list) != 0):
		#sorted_bed_list.sort(key=operator.attrgetter('start'));
		total_number_tags = len(sorted_bed_list);
		current_tag = sorted_bed_list[0];
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count>=threshold):
					current_tag.score = current_count;
					multiple_copy_read_list.append(current_tag);
				current_tag = item;
				current_value = item.start;
				current_count = 1; #reset
			else:
				current_count += 1;
		#last read	
		if (current_count>=threshold):
			item.score = current_count;
			multiple_copy_read_list.append(item);
	return multiple_copy_read_list;
	
	
	

def filter_reads(sorted_bed_list, cutoff, outfile):
	"""
	The histogram of tag copy will provide a cutoff for filtering the raw bed, as some of the tags has too many copies.  
	
	return filtered bed objects.
	"""
	if (len(sorted_bed_list) != 0):
		out = open(outfile, 'w')
		#sorted_bed_list.sort(key=operator.attrgetter('start'));
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_end = (sorted_bed_list[0]).end;
		current_count = 1;
		current_tag = sorted_bed_list[0];
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count <= cutoff):
					write(current_tag, out);
				current_value = item.start;
				current_end = item.end;
				current_count = 1; 
				current_tag = item;
			elif (item.end != current_end):
				if (current_count <= cutoff):
					write(current_tag, out);
				current_value = item.start;
				current_end = item.end;
				current_count = 1; 
				current_tag = item;
			else:
				if (current_count <= cutoff):
					write(current_tag, out);
				current_count += 1;
		if (current_count <= cutoff): #last tag
			write(current_tag, out);		
		out.close();


def filter_reads_add(sorted_bed_list, cutoff, outfile):
	"""
	The histogram of tag copy will provide a cutoff for filtering the raw bed, as some of the tags has too many copies.  
	
	return filtered bed objects.
	"""
	if (len(sorted_bed_list) != 0):
		out = open(outfile, 'a')
		#sorted_bed_list.sort(key=operator.attrgetter('start'));
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_end = (sorted_bed_list[0]).end;
		current_count = 1;
		current_tag = sorted_bed_list[0];
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count <= cutoff):
					write(current_tag, out);
				current_value = item.start;
				current_end = item.end;
				current_count = 1; 
				current_tag = item;
			elif (item.end != current_end):
				if (current_count <= cutoff):
					write(current_tag, out);
				current_value = item.start;
				current_end = item.end;
				current_count = 1; 
				current_tag = item;
			else:
				if (current_count <= cutoff):
					write(current_tag, out);
				current_count += 1;
				
		if (current_count <= cutoff): #last tag
			write(current_tag, out);		
		out.close();


def write (item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + item.name + "\t" + str(int(item.score)) + "\t" + item.strand + "\n";
	out.write(outline);

	
def write_list (bed6_list, out):
	"""
	write a bed6_list into outfile. The file openning and closing is handled from outside. 
	"""
	for item in bed6_list:
		#chrom, start, end, name, score, strand
		outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + item.name + "\t" + str(int(item.score)) + "\t" + item.strand + "\n";
		out.write(outline);
	
	
def combine_histogram(a, b):
	t=[];
	if len(a)<len(b):
		t = b;
		for index in xrange(len(a)):
			t[index]  += a[index];
	else:
		t = a;
		for index in xrange(len(b)):
			t[index]  += b[index];
	return t;

def write_histogram(a, outfile):
	out = open(outfile, 'w');
	for index in xrange(len(a)):
		if (a[index] != 0):
			outline = str(index) + "\t" + str(a[index]) +"\n";
			out.write(outline); 
	out.close();
	
def combine_read_copy_distribution(species, file_name):
	"""
	file_name is for the raw tag file. 
	need BED6 to split the positive and negative tags. 
	"""
	chroms = GenomeData.species_chroms[species];
	histogram =[];
	bed_vals = BED.BED(species, file_name, "BED6", 0);
	for chrom in chroms:
		if chrom in bed_vals.keys():
			sorted_bed_list = (bed_vals[chrom]).sort(key=operator.attrgetter('start'));
			(plus_bed_list, minus_bed_list) = breakUpStrands(sorted_bed_list);
			plus_histogram = find_read_copy_distribution(plus_bed_list);
			histogram = combine_histogram(plus_histogram, histogram)
			minus_histogram = find_read_copy_distribution(minus_bed_list);
			histogram = combine_histogram(minus_histogram, histogram);
	return histogram;	

def find_reads_in_region(start, end, reads):
	"""
	Find all the raw mapped reads in a given interval
	reads: a list of bed6 objects.
	Return the the list of reads in the region. 
	There is an alternative method written by Chongzhi, which perhaps is better.
	"""
	reads.sort(key=operator.attrgetter('start'));
	read_starts = [];
	for item in reads:
		read_starts.append(item.start);
	
	# if start is not in the read_starts, a[i}<s<a[i+1], then start_ind will be i+1, ok, 
	# if start is in the read_starts, s = a[i]=a[i+1]=,,,, start_ind = i, ok
	start_ind = bisect.bisect_left[start, read_starts];
	# if end is not in the read_starts, a[i]<e<a[i+1], then end_ind = i, ok
	# if end is in the read_starts, a[i-1]= a[i] = end, end_ind = i 
	end_ind = bisect.bisect_right[end, read_starts]-1;
	return reads[start_ind: end_ind+1];

def separateByChrom(chroms, file):
    for chrom in chroms:
        match = chrom + "[[:space:]]";
        tmpFile = chrom + ".bed";
        try:
            if os.system('%s %s %s > %s' %
                         (grep, match, file, tmpFile)): raise
        except: sys.stderr.write("grep failed\n");


def tag_position(sline, fragment_size):
	shift = int(round(fragment_size/2))
	if plus.match(sline[5]):
		return atoi(sline[1]) + shift
	elif minus.match(sline[5]):
		return atoi(sline[2]) - shift
	
	
def filter_tags_by_islands(chroms, islands, bed_file, filtered_file, fragment_size):
	"""
	Input:
		chroms a list of chrom
		islands is a dictionary keyed by chrom. 
		bed_file is for the raw reads.
			whether a tag belongs to an islands depends on the fragment_size
		filtered_file is the output file for the filtered reads.
		fragment_size results in the raw reads shift.
	Return 
		list of tags on the islands 
		list of tags not on the islands
	"""
	for chrom in chroms:
		if chrom in islands.keys():
			island_start_list = []
			island_end_list = []
			for item in islands[chrom]:
				island_start_list.append(item.start)
				island_end_list.append(item.end)
			f = open(bed_file,'r')
			o = open(filtered_file, 'w')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					position = tag_position(sline, fragment_size)
					if bisect.bisect_right(island_start_list, position) - bisect.bisect_left(island_end_list, position) == 1:
						o.write('\t'.join(sline)+'\n')
			f.close()
			o.close()

			
			
def combineAllGraphFiles(chroms, final_out):
	"""
		Combine the seperately processed chromosomes, return the output file name
	"""
    	outfile = open(final_out,'w');
    	outfile.close();
    
    	for chrom in chroms:
            graph_file = chrom + "_filtered.bed";
            try:
                if os.system('%s %s >> %s' %
                         	(cat, graph_file, final_out)): raise
            except: sys.stderr.write("cat failed\n")
	return final_out
		

def cleanup(chroms):
    for chrom in chroms:
        bed_file = chrom + ".bed";
        graph_file = chrom + "_filtered.bed";
        try:
            if os.remove('%s' %
                         (bed_file)): raise
        except: sys.stderr.write("clean up failed\n");
        try:
            if os.remove('%s' %
                         (graph_file)): raise
        except: sys.stderr.write("clean up failed\n");

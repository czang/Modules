#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import separate_by_chrom


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This function is used in island comparison of different cell types. The plan here is to decide whether the overlapping island is similar. The two islands are already determined to be overlapping by other function. The similarity is determined by two criterions
1) the tag counts in the overlapped area should be at least x percent of the total count of the respective islands.
2) the ratio of the normalized tag counts in the overlapped area between the two islands has to be close to 1.

'''


def getTagCounts(bed, chrom):
        '''
		bed should be a bed object with the bed_graph data type
		The data structure is a dictionary keyed by chromosome and valued by   
				a bunch of bed_graph object each of which contains [chrom start end count]
        Return a list of TagCounts on a given chromosome
        '''
        TagCounts = [];
        for t in bed.bed_vals[chrom]:
            TagCounts.append(t.value);
        try:
            return TagCounts;
        except:
            sys.stderr.write("Having trouble returning ends %s\n" % TagCounts)
            return ''


def FindTagCountInOverlap (OverlapRegion, StartList, EndList, TagCountList):
        ''' Return the tag counts in the overlapped region
                OverlapRegion is a dictionary
                The others are lists

                Returns the tag count in the overlapped region
                
                StartPosition is where the OverlapRegion[start] is supposed to be. 
                Since OverlapRegion[start] is the start of  the overlap region, 
				it has to be equal or bigger than the start of the first island, 
				and smaller or equal to the end of the last island. But it could be
				bigger than the start of the last island! In this case, the returned 
				value of the indexsearch is len(list), which is out of the range
		of the list.
                
        '''
        # StartPosition is where the OverlapRegion[start] is supposed to be
	StartPosition = bisect.bisect_left(StartList, OverlapRegion['start'])
	if StartPosition <len(StartList):
		#assert OverlapRegion['start'] <= StartList[StartPosition] and OverlapRegion['start']>StartList[StartPosition-1]
		if StartPosition == 0:
			i = 0;
        	elif  OverlapRegion['start'] <= EndList[StartPosition-1]: # The first bin is the one with index StartPosition-1
        		i = StartPosition-1;
        	else:  #The first bin is the one with index StartPosition
                	i = StartPosition;
           	
		TotalCount = 0.0;
        	while OverlapRegion['end'] >= StartList[i] :
        		TotalCount += TagCountList[i];
            		i = i+1;
			if i == len(StartList): break;

	elif StartPosition == len(StartList):
		TotalCount = TagCountList[-1];
		
	return TotalCount;


def aver_tag_density_on_islands(islands):
	'''return the average tag density of the islands in one chromosome. It is the ratio between the total tag number in all islands and the total island lengths '''
	total_bp = 0.0
	total_counts = 0.0
	for island_index in range(0, len(islands)):
		total_bp = total_bp - islands[island_index].start + islands[island_index].end + 1
		total_counts = total_counts + islands[island_index].value
	return float(total_counts) / float(total_bp)


def find_islands_overlapping(bed_vals_1, bed_vals_2, chrom, f, g, h, k):
	if chrom in bed_vals_1.keys() and chrom in bed_vals_2.keys():
		
		islands_1 = bed_vals_1[chrom]
		islands_2 = bed_vals_2[chrom]
		m = 0 # Number of real overlapped pairs.
		#n = 0
		#p = 0 # Number of overlapped pairs, before determination with any creteria
		nh = 0 # Number of NOT overlapped islands in islands_1
		nk = 0 # Number of NOT overlapped islands in islands_2
		Overlap_islands_numbers = {}
		if len(islands_1) != 0 and len(islands_2) != 0:
			islands_2_start = []
			islands_2_end = []
			for j in range(0, len(islands_2)):
				islands_2_start.append(islands_2[j].start)
				islands_2_end.append(islands_2[j].end)
			OverlapRegion = {}
			overlap_islands_2_start = []
			for island_index in range(0, len(islands_1)):
				index_start = bisect.bisect_left(islands_2_end, islands_1[island_index].start)
				index_end = bisect.bisect_left(islands_2_start, islands_1[island_index].end)
				''' Tricky point: search end in start list and search start in end list '''
				if index_end > index_start:  # Overlap determination
					#n+=1; # Number of overlapped islands in islands_1, before determination with any creteria. 
					z = 0
					for i in range(index_start, min(index_end, len(islands_2_end))):
						#p+=1  # Number of overlapped pairs, before determination with any creteria
						f.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) + '\n')
						g.write(chrom + ' ' + str(islands_2[i].start) + ' ' + str(islands_2[i].end) + '\n')
						m += 1 # Number of real overlapped pairs. 
						z += 1
						overlap_islands_2_start.append(islands_2[i].start)
					if z == 0:
						h.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) + '\n')
						nh+=1
				else:
					h.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) +'\n')
					nh+=1 
		
			for i in range (0, len(islands_2_start)): 
				if not (islands_2_start[i] in overlap_islands_2_start):
					k.write(chrom + ' ' + str(islands_2[i].start) + ' ' + str(islands_2[i].end) + '\n')
					nk+=1
			Overlap_islands_numbers['overlap1'] = len(islands_1) - nh
			Overlap_islands_numbers['overlap2'] = len(islands_2) - nk
			Overlap_islands_numbers['islands1'] = len(islands_1)
			Overlap_islands_numbers['islands2'] = len(islands_2)
			Overlap_islands_numbers['pairs'] = m  #pairs
			print chrom, 'total number:',len(islands_1),',',len(islands_2),'. Overlaped pairs:', m, ',', float(Overlap_islands_numbers['overlap1'])/float(len(islands_1))*100,'%,',float(Overlap_islands_numbers['overlap2'])/float(len(islands_2))*100,'%.'##'Average densities:', aver_tag_density_1, ',', aver_tag_density_2
		else:
			Overlap_islands_numbers['overlap1'] = 0
			Overlap_islands_numbers['overlap2'] = 0
			Overlap_islands_numbers['islands1'] = 0
			Overlap_islands_numbers['islands2'] = 0
			Overlap_islands_numbers['pairs'] = 0  #pairs
			print chrom, 'total number:',len(islands_1),',',len(islands_2),'. No Overlaped islands.'
	return Overlap_islands_numbers

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--bedfile1", action="store", type="string", dest="bedfile1", metavar="<file>", help="file 1 with islands info to be compared")
	parser.add_option("-b", "--bedfile2", action="store", type="string", dest="bedfile2", metavar="<file>", help="file 2 with islands info to be compared")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	parser.add_option("-p", "--result1", action="store", type="string", dest="result1", metavar="<file>", help="result file name for bed file 1")
	parser.add_option("-q", "--result2", action="store", type="string", dest="result2", help="result file name for bed file 2", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	total_overlap_number_1 = 0
	total_overlap_number_2 = 0
	total_islands_1 = 0
	total_islands_2 = 0
	total_pairs = 0
	
	separate_by_chrom.separateByChrom(chroms, opt.bedfile1, '.island1')
	separate_by_chrom.separateByChrom(chroms, opt.bedfile2, '.island2')
	
	
	for chrom in chroms: 
		f = open(chrom + '.similar1', 'w')
		g = open(chrom + '.similar2', 'w')
		h = open(chrom + '.different1', 'w')
		k = open(chrom + '.different2', 'w')
		bed_vals_1 = BED.BED(opt.species, chrom+'.island1', "BED3", 0)
		bed_vals_2 = BED.BED(opt.species, chrom+'.island2', "BED3", 0)
		overlap_number = find_islands_overlapping(bed_vals_1, bed_vals_2, chrom, f, g, h, k)
		total_overlap_number_1 = total_overlap_number_1 + overlap_number['overlap1']
		total_overlap_number_2 = total_overlap_number_2 + overlap_number['overlap2']
		total_islands_1 = total_islands_1 + overlap_number['islands1']
		total_islands_2 = total_islands_2 + overlap_number['islands2']
		total_pairs = total_pairs + overlap_number['pairs']
		f.close()
		g.close()
		h.close()
		k.close()
	
	print total_overlap_number_1,'/',total_islands_1, float(total_overlap_number_1)/float(total_islands_1)*100, '%;', total_overlap_number_2, '/',total_islands_2,float(total_overlap_number_2)/float(total_islands_2)*100,'%.'
	
	separate_by_chrom.combineAllGraphFiles(chroms, '.similar1', opt.result1+'_overlap.island')
	separate_by_chrom.combineAllGraphFiles(chroms, '.similar2', opt.result2+'_overlap.island')
	separate_by_chrom.combineAllGraphFiles(chroms, '.different1', opt.result1+'_nonoverlap.island')
	separate_by_chrom.combineAllGraphFiles(chroms, '.different2', opt.result2+'_nonoverlap.island')
	
	separate_by_chrom.cleanup(chroms, '.similar1')
	separate_by_chrom.cleanup(chroms, '.similar2')
	separate_by_chrom.cleanup(chroms, '.different1')
	separate_by_chrom.cleanup(chroms, '.different2')
	separate_by_chrom.cleanup(chroms, '.island1')
	separate_by_chrom.cleanup(chroms, '.island2')

if __name__ == "__main__":
	main(sys.argv)
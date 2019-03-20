#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

# get BED module
import BED

#Changes made by Weiqun: output files go to the current directory, values are interger.

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

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


def indexsearch(list, search):    
    """ classic binary search -- If 'search' is in 'list', its index is
    returned.  If not, the index where 'search' should be in 'list' is returned """
    right = len(list)
    left = 0
    previous_center = 0
    if search < list[0]:
        return 0
    while 1:
        center = (left + right) / 2
        candidate = list[center]
        if search == candidate:
            return center
        if center == previous_center:
            return (1 + center);
        elif search < candidate:
            right = center
        else:
            left = center
        previous_center = center


'''
This function is used in island comparison of different cell types. The plan here is to decide whether the overlapping island is similar. The two islands are already determined to be overlapping by other function. The similarity is determined by two criterions
1) the tag counts in the overlapped area should be at least x percent of the total count of the respective islands.
2) the ratio of the normalized tag counts in the overlapped area between the two islands has to be close to 1.

'''

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
	StartPosition = indexsearch(StartList, OverlapRegion['start'])
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


def find_islands_overlapping(bed_vals_1, bed_vals_2, Histogram1_bed, Histogram2_bed, OverlapThreshold, Significance, chrom, f, g, h, k):
	if chrom in bed_vals_1.keys() and chrom in bed_vals_2.keys() and chrom in Histogram1_bed.keys() and chrom in Histogram2_bed.keys():
		StartList1 = Histogram1_bed.getStarts(chrom)
		EndList1 = Histogram1_bed.getEnds(chrom)
		StartList2 = Histogram2_bed.getStarts(chrom)
		EndList2 = Histogram2_bed.getEnds(chrom)
		TagCountList1 = getTagCounts(Histogram1_bed, chrom)
		TagCountList2 = getTagCounts(Histogram2_bed, chrom)
		
		islands_1 = bed_vals_1[chrom]
		islands_2 = bed_vals_2[chrom]
		m = 0 # Number of real overlapped pairs.
		#n = 0
		#p = 0 # Number of overlapped pairs, before determination with any creteria
		nh = 0 # Number of NOT overlapped islands in islands_1
		nk = 0 # Number of NOT overlapped islands in islands_2
		Overlap_islands_numbers = {}
		if len(islands_1) != 0 and len(islands_2) != 0:
			aver_tag_density_1 = aver_tag_density_on_islands(islands_1)
			aver_tag_density_2 = aver_tag_density_on_islands(islands_2)
			
			islands_2_start = []
			islands_2_end = []
			for j in range(0, len(islands_2)):
				islands_2_start.append(islands_2[j].start)
				islands_2_end.append(islands_2[j].end)
			OverlapRegion = {}
			overlap_islands_2_start = []
			for island_index in range(0, len(islands_1)):
				index_start = indexsearch(islands_2_end, islands_1[island_index].start)
				index_end = indexsearch(islands_2_start, islands_1[island_index].end)
				''' Tricky point: search end in start list and search start in end list '''
				if index_end > index_start:  # Overlap determination
					#n+=1; # Number of overlapped islands in islands_1, before determination with any creteria. 
					z = 0
					for i in range(index_start, min(index_end, len(islands_2_end))):
						#p+=1  # Number of overlapped pairs, before determination with any creteria
						OverlapRegion['start'] = max (islands_1[island_index].start , islands_2[i].start)
						OverlapRegion['end'] = min (islands_1[island_index].end , islands_2[i].end)
						Overlap_Length = float(OverlapRegion['end'] - OverlapRegion['start'] +1 )
						OverlapCounts_1 = FindTagCountInOverlap (OverlapRegion, StartList1, EndList1, TagCountList1)
						OverlapConfirmation = -1
						#print OverlapCounts_1, OverlapRegion['start'], OverlapRegion['end']
						if OverlapCounts_1 >= float(islands_1[island_index].value) * OverlapThreshold :
							OverlapCounts_2 = FindTagCountInOverlap (OverlapRegion, StartList2, EndList2, TagCountList2)
							#print OverlapCounts_2
							if OverlapCounts_2 >= float (islands_2[i].value) * OverlapThreshold:
								Normalized_Counts_1 = float(OverlapCounts_1) / (Overlap_Length * aver_tag_density_1)
								Normalized_Counts_2 = float(OverlapCounts_2) / (Overlap_Length * aver_tag_density_2)
								ratio = abs(Normalized_Counts_1 - Normalized_Counts_2) / (Normalized_Counts_1 + Normalized_Counts_2)
								if ratio <= Significance:
									OverlapConfirmation = 1
						if OverlapConfirmation > 0:  # Real overlapping determination
							f.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) + ' ' + str(int(islands_1[island_index].value))+'\n')
							g.write(chrom + ' ' + str(islands_2[i].start) + ' ' + str(islands_2[i].end) + ' ' + str(int(islands_2[i].value)) + '\n')
							#print islands_1[island_index].start, '\t', islands_1[island_index].end, '\t', islands_1[island_index].value, '\t', Normalized_Counts_1, '\t', islands_2[i].start, '\t', islands_2[i].end, '\t', islands_2[i].value, '\t', Normalized_Counts_2
							m += 1 # Number of real overlapped pairs. 
							z += 1
							overlap_islands_2_start.append(islands_2[i].start)
					if z == 0:
						h.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) + ' ' + str(int(islands_1[island_index].value))+'\n')
						nh+=1
				else:
					h.write(chrom + ' ' + str(islands_1[island_index].start) + ' ' + str(islands_1[island_index].end) + ' ' + str(int(islands_1[island_index].value))+'\n')
					nh+=1 
		
			for i in range (0, len(islands_2_start)): 
				if not (islands_2_start[i] in overlap_islands_2_start):
					k.write(chrom + ' ' + str(islands_2[i].start) + ' ' + str(islands_2[i].end) + ' ' + str(int(islands_2[i].value)) + '\n')
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
	parser.add_option("-d", "--bedfile3", action="store", type="string", dest="Histogram1", metavar="<file>", help="file 1 with summary info in bed graph format")
	parser.add_option("-g", "--bedfile4", action="store", type="string", dest="Histogram2", metavar="<file>", help="file 2 with summary info in bed graph format")
	parser.add_option("-t", "--OverlapThreshold", action="store", type="float", dest="OverlapThreshold",  metavar="<float>", help="real number between 0 and 1")
	parser.add_option("-r", "--ratio", action="store", type="float", dest="Significance", metavar="<float>", help="real number between 0 and 1")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8 or hg18", metavar="<str>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
	
	mm_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX', 'chrY', 'chrM']
	hg_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
	
	if opt.species == "mm8": 
		chroms = mm_chroms
	elif opt.species == "hg18": 
		chroms = hg_chroms
	else:
		print 'Species Error'

	total_overlap_number_1 = 0
	total_overlap_number_2 = 0
	total_islands_1 = 0
	total_islands_2 = 0
	total_pairs = 0

	bed_vals_1 = BED.BED(opt.species, opt.bedfile1, "BED_GRAPH", 0)
	bed_vals_2 = BED.BED(opt.species, opt.bedfile2, "BED_GRAPH", 0)
	Histogram1_bed = BED.BED(opt.species, opt.Histogram1, "BED_GRAPH")
	Histogram2_bed = BED.BED(opt.species, opt.Histogram2, "BED_GRAPH")
	
	for chrom in chroms: 
		rr = (opt.bedfile1.split('/'))[-1]
		qq = (opt.bedfile2.split('/'))[-1]
		f = open(rr + '-' + str(chrom) + '-similar', 'w')
		g = open(qq + '-' + str(chrom) + '-similar', 'w')
		h = open(rr + '-' + str(chrom) + '-different', 'w')
		k = open(qq + '-' + str(chrom) + '-different', 'w')
		
		overlap_number = find_islands_overlapping(bed_vals_1, bed_vals_2, Histogram1_bed, Histogram2_bed, opt.OverlapThreshold, opt.Significance, chrom, f, g, h, k)
		total_overlap_number_1 = total_overlap_number_1 + overlap_number['overlap1']
		total_overlap_number_2 = total_overlap_number_2 + overlap_number['overlap2']
		total_islands_1 = total_islands_1 + overlap_number['islands1']
		total_islands_2 = total_islands_2 + overlap_number['islands2']
		total_pairs = total_pairs + overlap_number['pairs']
	
		f.close()
		g.close()
		h.close()
		k.close()
	#l = open(str(opt.bedfile1) + '-' + str(opt.bedfile2) + '-' + str(opt.OverlapThreshold) + '-' + str(opt.Significance),'w')
	l = open("Th1-Th2_islands" + str(opt.OverlapThreshold) + '-' + str(opt.Significance),'w')
	l.write(str(opt.OverlapThreshold) + '\t' + str(opt.Significance) + '\t' + str( total_overlap_number_1)+'\t'+str(total_islands_1)+'\t'+ str( float(total_overlap_number_1)/float(total_islands_1)) + '\t'+ str(total_overlap_number_2)+ '\t'+str(total_islands_2)+'\t' + str(float(total_overlap_number_2)/float(total_islands_2)))
	l.close()
	print 'OverlapThreshold:',opt.OverlapThreshold,', TagDifference:',opt.Significance,'.', total_overlap_number_1,'/',total_islands_1, float(total_overlap_number_1)/float(total_islands_1)*100, '%;', total_overlap_number_2, '/',total_islands_2,float(total_overlap_number_2)/float(total_islands_2)*100,'%.'


if __name__ == "__main__":
	main(sys.argv)
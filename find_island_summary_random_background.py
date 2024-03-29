#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
import GenomeData;
import associate_tags_with_regions
import SeparateByChrom
import get_total_tag_counts
import Utility
import scipy.stats

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--rawchipreadfile", action="store", type="string", dest="chipreadfile", metavar="<file>", help="raw read file from chip in bed format")
	parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>", help="average size of a fragment after CHIP experiment")
	parser.add_option("-d", "--islandfile", action="store", type="string", dest="islandfile", metavar="<file>", help="island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count summary file")
	parser.add_option("-t", "--mappable_faction_of_genome_size ", action="store", type="float", dest="fraction", help="mapable fraction of genome size", metavar="<float>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species];
		genomesize = sum (GenomeData.species_chrom_lengths[opt.species].values());
		genomesize = opt.fraction * genomesize;
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	chip_library_size=get_total_tag_counts.get_total_tag_counts(opt.chipreadfile);
	#control_library_size=get_total_tag_counts.get_total_tag_counts(opt.controlreadfile);
	print "chip library size  ", chip_library_size
	#print "control library size  ", control_library_size
	
	totalchip = 0;
	#totalcontrol = 0;
	
	islands = BED.BED(opt.species, opt.islandfile, "BED3", 0);
	
	# separate by chrom the chip library
	if Utility.fileExists(opt.chipreadfile):
		SeparateByChrom.separateByChrom(chroms, opt.chipreadfile, '.bed1');
	else:
		print opt.chipreadfile, " not found";
		sys.exit(1)
	# separate by chrom the control library
	#if Utility.fileExists(opt.controlreadfile):
	#	SeparateByChrom.separateByChrom(chroms, opt.controlreadfile, '.bed2');
	#else:
	#	print opt.controlreadfile, " not found";
	#	sys.exit(1)	
	
	island_chip_readcount = {};
	#island_control_readcount = {};
	
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				if Utility.is_bed_sorted(island_list) == 0:
					island_list.sort(key=operator.attrgetter('start'));
					
				island_start_list = []
				island_end_list = []
				for item in island_list:
					island_start_list.append(item.start)
					island_end_list.append(item.end)
	
				island_chip_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed1";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index =associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_chip_readcount_list[index] += 1;
							totalchip += 1;
				f.close();
				island_chip_readcount[chrom] = island_chip_readcount_list;
							
				'''island_control_readcount_list=[0]*len(island_list);
				read_file = chrom + ".bed2";
				f = open(read_file,'r')
				for line in f:
					if not re.match("#", line):
						line = line.strip()
						sline = line.split()
						position = associate_tags_with_regions.tag_position(sline, opt.fragment_size)
						index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position);
						if index >= 0:
							island_control_readcount_list[index] += 1;
							totalcontrol += 1;
				f.close();
							
				island_control_readcount[chrom] = island_control_readcount_list;'''			
						
	chip_background_read = chip_library_size - totalchip;
	#control_background_read = control_library_size - totalcontrol;
	#scaling_factor = chip_background_read*1.0/control_background_read;
	#scaling_factor = chip_library_size*1.0/control_library_size;
	
	
	print "Total number of chip reads on islands is: ", totalchip; 
	#print "Total number of control reads on islands is: ", totalcontrol; 

	#print "chip_background_read   ", chip_background_read
	#print "control_background_read   ", control_background_read

	out = open(opt.out_file, 'w');
	pvalue_list = [];
	result_list = [];
	for chrom in chroms:
		if chrom in islands.keys():
			if len(islands[chrom]) != 0:
				island_list = islands[chrom];
				for index in xrange(len(island_list)):
					item = island_list[index];
					length = item.end - item.start + 1;
					observation = (island_chip_readcount[chrom])[index];
					average = length * chip_library_size *1.0/genomesize;
					fc = float(observation)/float(average);
					'''#control_tag = (island_control_readcount[chrom])[index];
					if (island_control_readcount[chrom])[index] > 0:
						#average = (island_control_readcount[chrom])[index] * scaling_factor;
						average = control_tag * scaling_factor
						fc = float(observation)/float(average);
					else:
						length = item.end - item.start + 1;
						average = length * control_library_size *1.0/genomesize;			
						average = min(0.25, average)* scaling_factor;
						fc = float(observation)/float(average);'''
					if observation > average:
						pvalue =	 scipy.stats.poisson.sf(observation, average)[()]; 
					else:
						pvalue = 1.0;
					pvalue_list.append(pvalue);
					item_dic = {}
					item_dic['chrom'] = item.chrom
					item_dic['start'] = item.start
					item_dic['end'] = item.end
					item_dic['chip'] = observation
					#item_dic['control'] = control_tag
					item_dic['pvalue'] = pvalue
					item_dic['fc'] = fc
					result_list.append(item_dic)
				
	#pvalue_list.sort()
	for item in result_list:
		pvalue = float(item['pvalue'])
		#alpha = pvalue * len(result_list) / (pvalue_list.index(pvalue) + 1)
		outline = item['chrom'] + "\t" + str(item['start']) + "\t" + str(item['end']) + "\t" + str(item['chip']) + "\t0\t" + str(item['pvalue']) + "\t" + str(item['fc']) + "\n";	
		out.write(outline);		
	out.close();
	
	
	SeparateByChrom.cleanup(chroms, '.bed1');
	#SeparateByChrom.cleanup(chroms, '.bed2');
	
	

if __name__ == "__main__":
	main(sys.argv)

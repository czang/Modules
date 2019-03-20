#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED
import get_total_tag_counts
import Background_island_calculation
from species_chroms import *


""" 
	
	Take in coords for bed_gaph type summary files and find 'islands' of modifications.

	There are a number options here that can be turned on or off depending on need

		Right now:

	(1) scan for all 'islands' where the single window or consecutive windows

	(2) remove all single window 'islands' -- this can (should?) be commented
		out when looking for more localized signals (such as TF binding sites?)

	(3) try to combine islands that are within gap distance of each other.
		This gap distance is supplemented by a window_buffer just so we don't do
		anything stupid with window sizes

	(4) Remove all single window combined islands -- if step (2) was done,
		this is redundant

	(5) Lastly, filter out all the islands we've found that have a total
		score < islands_minimum_tags
"""

##  parameter
window_buffer = 10;

def combineProximalIslands(islands, gap):
    """
    Extend the islands found with the findIslands function
    """
    proximal_island_dist = gap + window_buffer;
    combined_islands = [];
    island_index = 0;
    got_last = 0;
  
    while island_index < (len(islands) - 1):
        i = islands[island_index];
        next_i = islands[island_index+1];
        start = i.start;
        end = i.end;
        score = i.value;
        extend = 1;
        while extend:
            if island_index < (len(islands) - 1):
                i = islands[island_index];
                next_i = islands[island_index+1];
                next_dist = abs(next_i.start - i.end);
                if next_dist <= proximal_island_dist:
                    end = next_i.end;
                    score += next_i.value;
                    island_index += 1;
                else:
                    extend = 0;
            else:
                extend = 0;
                got_last = 1;
        island_index += 1;
        whole_island = BED.BED_GRAPH(i.chrom, start, end, score);
        combined_islands.append(whole_island);
    if got_last == 0:
        combined_islands.append(islands[-1]);
    return combined_islands;


def removeSingleWindowIslands(islands, window):
    filtered_islands = [];
    for i in islands:
        size = i.end - i.start;
        if size > window:
            filtered_islands.append(i);
    return filtered_islands;


def findIslands(bed_vals, chrom, outfilename, window, gap, islands_minimum_tags):
    """
    Find all islands, which are _consecutive_ windows
    scoring about threshold
    """
    islands = [];
    if chrom in bed_vals.keys():
        tags = bed_vals[chrom];
        tags.sort(key=operator.attrgetter('start'));
        tag_index = 0;
        got_last = 0;
        while tag_index < (len(tags) - 1):
            island_start = tags[tag_index].start;
            island_end = tags[tag_index].end;
            island_value = tags[tag_index].value;
            extend = 1;
            while extend:
                if tag_index < (len(tags) - 1):
                    next_distance = abs(tags[tag_index+1].start - tags[tag_index].end);
                    if next_distance < window_buffer:
                        island_end = tags[tag_index+1].end
                        island_value += tags[tag_index+1].value
                        tag_index += 1;
                    else:
                        island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
                        islands.append(island);
                        extend = 0;
                else:
                    got_last = 1;
                    island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
                    islands.append(island);
                    extend = 0;
            tag_index += 1;

        """ add last on if didn't get there """
        if got_last == 0:
            island_start = tags[tag_index].start;
            island_end = tags[tag_index].end;
            island_value = tags[tag_index].value;

            island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
            islands.append(island);


        """ remove all single window islands """
        #islands = removeSingleWindowIslands(islands, window);

        """ combine all islands within """
        combined_islands = combineProximalIslands(islands, gap);

        """ remove all single window islands """
        combined_islands = removeSingleWindowIslands(combined_islands, window);

        
        """ filter islands by tag counts and print out """
        outfile = open(outfilename, 'w');
        for i in combined_islands:
            if i.value >= islands_minimum_tags:
                outline = chrom + " " + str(i.start) + " " + str(i.end) + " " + str(i.value) + "\n";
                outfile.write(outline);
        outfile.close();


def combine_proximal_regions(islands, gap, window_size_buffer=10):
    """
    	Extend the regions found in the find_continuous_region function
    	if gap is not allowed, gap = 0, if one window is allowed, gap = window_size (200)
	return a list of combined regions. 
    """
    proximal_island_dist = gap + window_size_buffer;
    combined_islands = [];
    island_index = 0;
    got_last = 0;
  
    while island_index < (len(islands) - 1):
        i = islands[island_index];
        next_i = islands[island_index+1];
        start = i.start;
        end = i.end;
        score = i.value;
        extend = 1;
        while extend:
            if island_index < (len(islands) - 1):
                i = islands[island_index];
                next_i = islands[island_index+1];
                next_dist = abs(next_i.start - i.end);
                if next_dist <= proximal_island_dist:
                    end = next_i.end;
                    score += next_i.value;
                    island_index += 1;
                else:
                    extend = 0;
            else:
                extend = 0;
                got_last = 1;
        island_index += 1;
        whole_island = BED.BED_GRAPH(i.chrom, start, end, score);
        combined_islands.append(whole_island);
    if got_last == 0:
        combined_islands.append(islands[-1]);
    return combined_islands;


def find_region(bed_vals, chrom, gap, window_size, window_size_buffer=10):
	"""
		bed_graph_file is the bed_graph type summary file
                
    		Find all regions  made of consecutive windows. The requirement for each window to score
                above tag_count_threshold_in_window is implemented by filtering involved in the making of bed_vals
    	
    		return a list of BED.BED_GRAPH, each of which has (i.chrom, start, end, score)
    	"""  
     
	islands = [];
	if chrom in bed_vals.keys():
		tags = bed_vals[chrom];
		if (len(tags) >0):
			tags.sort(key=operator.attrgetter('start'));
			tag_index = 0;
			got_last = 0;
			while tag_index < (len(tags) - 1):
				island_start = tags[tag_index].start;
				island_end = tags[tag_index].end;
				island_value = tags[tag_index].value;
				extend = 1;
				while extend:
					if tag_index < (len(tags) - 1):
						next_distance = abs(tags[tag_index+1].start - tags[tag_index].end);
						if next_distance < window_size_buffer:
							island_end = tags[tag_index+1].end
							island_value += tags[tag_index+1].value
							tag_index += 1;
						else:
							island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
							islands.append(island);
							extend = 0;
					else:
						got_last = 1;
						island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
						islands.append(island);
						extend = 0;
				tag_index += 1;

			# add last on if didn't get there
			if got_last == 0:
				#print chrom, tag_index;
				island_start = tags[tag_index].start;
				island_end = tags[tag_index].end;
				island_value = tags[tag_index].value;
	
				island = BED.BED_GRAPH(chrom, island_start, island_end, island_value);
				islands.append(island);
	
			#combine all islands within 
			combined_islands =combine_proximal_regions(islands, gap, window_size_buffer);
			
			#""" remove all single window islands """
			#combined_islands = removeSingleWindowIslands(combined_islands, window_size);

			
		else: combined_islands = [];	
		return combined_islands;
	else:
		print "Chromosome number not right!!";


def find_region_above_threshold(island_list, islands_minimum_tags):
	filtered_islands = [];
	for island in island_list:
		if island.value >= islands_minimum_tags: filtered_islands.append(island);
	return filtered_islands;

def find_region_above_threshold_from_file(Islands_file, species, islands_minimum_tags, out_islands_file):
	if species_chroms.has_key(species):
		bed_vals = BED.BED(species, Islands_file, "BED_GRAPH");
		outputfile = open(out_islands_file, 'w');
		total_number_islands = 0.0;
		total_tags_on_islands = 0.0;
		
		for chrom in bed_vals.keys():
			islands = bed_vals[chrom];
			islands = find_region_above_threshold(islands, islands_minimum_tags);
			total_number_islands += len(islands);
			for i in islands:
                    		outline = chrom + " " + str(i.start) + " " + str(i.end) + " " + str(i.value) + "\n";	
                    		outputfile.write(outline);
				total_tags_on_islands += i.value;
            	outputfile.close();
		print "Total tag counts on filtered islands: ", total_tags_on_islands;
	else:
		print opt.species + " is not in the species list";


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
			dest="species", help="mm8, hg18, background", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bedfile",  metavar="<file>", 
			help="file with summary info in bed graph format")
	parser.add_option("-o", "--outfile", action="store", type="string",
			dest="outfile", help="output island file name", metavar="<file>")
	parser.add_option("-w", "--window_size", action="store", type="int",
			dest="window_size", help="window size of summary", metavar="<int>")
	parser.add_option("-g", "--gap_size", action="store", type="int",
			dest="gap", help="gap size (in bps) to allow", metavar="<int>")
    	parser.add_option("-m", "-- min_tag_count_in_window", action="store", type="int",
			dest="tag_count_threshold_in_window", help="minimum number of tags required in each window", metavar="<int>")
	parser.add_option("-t", "--min_tag_count_in_island", action="store", type="int",
			dest="islands_minimum_tags", help="minimum number of tags required in each island; if <0, self-determine according to background and p value; if>=0 use it", metavar="<int>")
	parser.add_option("-p", "--p_value_threshold", action="store", type="float",
			dest="p_val", help="p value threshold, .1", metavar="<float>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
		parser.print_help()
		sys.exit(1)

        if species_chroms.has_key(opt.species):
		total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.bedfile);
		genome_length = 0.0;
		for item in species_chrom_lengths[opt.species]: genome_length += item;
		average = opt.window_size * total_tag_counts/genome_length;
		scaled_genome_length = int(round(genome_length/opt.window_size));
		scaled_gap = int(round(opt.gap/opt.window_size));
		
		#print total_tag_counts, genome_length, scaled_genome_length, scaled_gap;
		
		if (opt.islands_minimum_tags<0):
			islands_minimum_tags = Background_island_calculation.find_island_threshold(opt.p_val, opt.tag_count_threshold_in_window, scaled_gap, average, scaled_genome_length);
		else: 	
			islands_minimum_tags = opt.islands_minimum_tags;


            	"""  Get all bed vals with a score >= window_count_threshold"""
            	bed_vals = BED.BED(opt.species, opt.bedfile, "BED_GRAPH", opt.tag_count_threshold_in_window);
		
	    	
		#islands_minimum_tags = 0;
	    	print opt.bedfile;
		print "Total tag counts: "+ str(total_tag_counts);
		print "Genome Length is: ", genome_length;
		print "Island threshold is " + str(islands_minimum_tags);
		
            	total_number_islands = 0.0;
		total_tags_on_islands = 0.0;
		outputfile = open(opt.outfile, 'w');
            	for chrom in bed_vals.keys():
			islands = find_region(bed_vals, chrom, opt.gap, opt.window_size, 10);
			islands = find_region_above_threshold(islands, islands_minimum_tags);
			total_number_islands += len(islands);
			
			for i in islands:
                    		outline = chrom + " " + str(i.start) + " " + str(i.end) + " " + str(i.value) + "\n";	
                    		outputfile.write(outline);
				total_tags_on_islands += i.value;
            	outputfile.close();

            	print "The total number of islands is " +str(total_number_islands);
		print "Coverage is ",  total_tags_on_islands/float(total_tag_counts);
        else:
            	print opt.species + " is not in the species list";

if __name__ == "__main__":
    	main(sys.argv)


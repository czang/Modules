#!/usr/bin/env python
import re, os, sys, shutil

bg_number_chroms = 10;
bg_length_of_chrom = 200000000;
background_chroms = [];
background_chrom_lengths = [];
for i in range(0,bg_number_chroms):
	background_chroms.append( 'chr' + str(i+1));
	background_chrom_lengths.append(bg_length_of_chrom);
#print background_chroms;



mm_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX', 'chrY', 'chrM']
hg_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
sacCer_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9', 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chrM'];
dm_chroms = ['chr2h', 'chr2L', 'chr2R', 'chr3h', 'chr3L', 'chr3R', 'chr4', 'chr4h', 'chrM', 'chrU', 'chrX', 'chrXh', 'chrYh'];


mm_chrom_lengths = [197069962,181976762, 159872112, 155029701, 152003063, \
	149525685, 145134094, 132085098, 124000669, 129959148, 121798632, \
	120463159, 120614378, 123978870, 103492577, 98252459, 95177420, 90736837, \
	61321190, 165556469, 16029404, 16299];
hg_chrom_lengths = [247249719, 242951149, 199501827, 191273063, 180857866, \
                    170899992, 158821424, 146274826, 140273252, 135374737, \
                    134452384, 132349534, 114142980, 106368585, 100338915, \
                    88827254, 78774742, 76117153, 63811651, 62435964, \
                    46944323, 49691432, 154913754, 57772954, 16571]  


species_chroms = {'mm8':mm_chroms, 'hg18':hg_chroms, "dm2":dm_chroms, "sacCer1":sacCer_chroms, 'background':background_chroms};
species_chrom_lengths={'mm8':mm_chrom_lengths, 'hg18':hg_chrom_lengths, 'background':background_chrom_lengths};


def main(argv):
    print species_chroms.keys();
    print;
    for species in species_chroms.keys():
        print species, ": ", species_chroms[species];
        print;

if __name__ == "__main__":
	main(sys.argv)


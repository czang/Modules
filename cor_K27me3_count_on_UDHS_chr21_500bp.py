#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import time
#import tables
import numpy as np
#from sklearn import cluster as sklearn_cluster
import h5py
import argparse,sys
from gene_set_manipulation import *


def read_hdf5( file_name, sample_names ):
    """ Apply motif stat function to all data in motif_file_name. 
    Data in numpy array val corresponds to idx entries. If idx if None all entries are used.""" 
    h5file = h5py.File(file_name)
    # X = None
    chr21_index=range(798590,808935)
    X = np.zeros( ( len(sample_names),len(chr21_index)) ) 
    for index,elem in enumerate(sample_names):
        a = h5file[elem]
        m = a.value
        X[index,:] = m[chr21_index]
        print index+1,"samples have been read"
    X = X.transpose()
    h5file.close()
    return X




def getSampleNames_hdf5(rp_hdf5):
    h5file = h5py.File(rp_hdf5)
    samplenames = []
    for sample in h5file.keys():
        if sample not in samplenames and sample not in ["chr","mid"]:
            samplenames.append(sample)
        else:
            continue
    h5file.close()
    return samplenames


def is_list_sorted(list):
        """
        Check if sorted in ascending order.
        input is a list of values.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(list)-1):
                if list[index] > list[index+1]:
                        sorted = 0;
        return sorted;


def getLocList(infile):
	f = open(infile,'r')
	List = []
	for line in f:
		line = line.strip()
		sline = line.split()
		start = atoi(sline[1])
		end = atoi(sline[2])
		loc = (start + end) / 2
		List.append(loc)
	f.close()
	assert is_list_sorted(List)
	return List


def cal_neighboring_cor(Data, corList, max_distance=500000):
	assert len(corList) == len(Data[:,0])
	x = []
	y = []
	z = []
	k = 0
	for i in range(0,len(corList)-2):
		j = i + 1
		while corList[j]-corList[i] <= max_distance:
			x.append(corList[i])
			y.append(corList[j])
			d = np.vstack((np.sqrt(Data[i]),np.sqrt(Data[j])))
			c = np.corrcoef(d)
			cor = c[0,1]
			z.append(round(cor,6))
			k += 1
			#print k
			if j < len(corList)-1:
				j += 1
			else:
				break
	return x, y, z


#hdf5_file="/data5/home/chenfei/JingyuFan/data_collection/MARGE/union_DHS_2016/hg38_H3K27ac_UDHS_1kb_count.h5"
#hdf5_file="/data5/home/chenfei/JingyuFan/data_collection/MARGE/union_DHS_2016/hg38_H3K27ac_1kb_count_UDHS_500.h5"
hdf5_file="/data5/home/chenfei/JingyuFan/data_collection/MARGE/union_DHS_2016/hg38_H3K27me3_1kb_count_UDHS_500.h5"

all_sample_name = getSampleNames_hdf5(hdf5_file)

with open("/data5/home/chenfei/JingyuFan/data_collection/MARGE/union_DHS_2016/hg38_H3K27me3_datasets.txt","r") as infile:
    temp=[]
    for line in infile:
        if line.strip():
            if line.strip() in all_sample_name:
                temp.append(line.strip())


all_sample_name = temp  
all_sample_name = get_gene_list("/data5/home/chenfei/JingyuFan/data_collection/MARGE/union_DHS_2016/hg38_H3K27me3_datasets.txt", 1)

start = time.time()
# all_sample_name = all_sample_name[:30]
data = read_hdf5(hdf5_file, all_sample_name)
#np.savetxt('test.txt', data)


DHS_file = "/data/home/czang/cistrome_integration/20160328_chr21_cormap/UDHS_chr21.bed"

loc_list =  getLocList(DHS_file)

(x, y, z) = cal_neighboring_cor(data, loc_list)

o = open("/data/home/czang/cistrome_integration/20160331_chr21_K27me3/chr21_map.txt",'w')
for i in range(0, len(x)):
	o.write(str(x[i]) + '\t' + str(y[i]) + '\t' + str(z[i]) + '\n')
o.close()

#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from gene_set_manipulation import *

cutoff = 15

genelist = get_gene_list(sys.argv[1]+'.txt',0)
f = open(sys.argv[1]+'_list.txt','w')
for item in genelist:
	item = item.strip()
	item = item.split('_')
	f.write(item[0]+'\n')
f.close()

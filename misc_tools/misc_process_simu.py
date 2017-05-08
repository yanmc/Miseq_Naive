#!/usr/bin/env python
# encoding: utf-8
"""
miss_trans_nucltide_on_hg19.py

Created by Mingchen on 2014-12-11.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  original fasta file
orgnism:  human,mouse or rabbit
referecce_file:  the file contain variable region start and end info eg: human_get_vdj.txt 

note: already trans the reversed reads 
"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import numpy as np
import matplotlib.pyplot as plt
from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
from misc_prepare_pbs import *
from misc_get_trimmed_region import *
from itertools import groupby
try:
    import cPickle as pickle
except ImportError:
    import pickle
from collections import Counter
def get_simu_dict(length, depth):
	repeat_region = csv.reader(open("/zzh_gpfs02/zhangyanfang/repeat_region_analysis/all_repeat_region_interval.txt", "rU"), delimiter = "\t")
	repeat_region_dict = {}
	for index, line in enumerate(repeat_region):
		repeat_region_dict.setdefault(line[0], []).append((line[1], line[2]))


	simu_record = csv.reader(open("/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record.txt"%(length, depth), "rU"), delimiter = "\t")

	simu_dict, repeatregion_simu_dict, non_repeatregion_simu_dict = {},{},{}
	print IUB_CODE
	for index, line in enumerate(simu_record):
		if index % 10000 == 0:
			print "Processed 10000 lines in simu %s %s."%(length, depth)
		#print line
		chromsome = line[0]
		position  = int(line[1])
		ref_nucle = line[2]
		snv_nucle = IUB_CODE[line[3]]
		haplo_type = line[4]
		in_repeat_region = 0
		try:
			for region in repeat_region_dict[chromsome]:
				if int(region[0]) <= position-1 <= int(region[1]):
					in_repeat_region = 1
					#print chromsome, position, ref_nucle, snv_nucle, haplo_type, in_repeat_region
					#sys.exit(0)
					break
		except KeyError:
			pass

		#print ref_nucle, snv_nucle, haplo_type
		simu_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
		if in_repeat_region == 1:
			repeatregion_simu_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
		else:
			non_repeatregion_simu_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
		#print simu_dict
		#if index == 10:
		#	sys.exit(0)

	#How many position be detected
	#print len(set(vcf_dict.keys())), len(set(simu_dict.keys())), len(set(vcf_dict.keys()) & set(simu_dict.keys()))
	#print len(set(repeatregion_vcf_dict.keys())), len(set(repeatregion_simu_dict.keys())), len(set(repeatregion_vcf_dict.keys()) & set(repeatregion_simu_dict.keys()))
	#print len(set(non_repeatregion_vcf_dict.keys())), len(set(non_repeatregion_simu_dict.keys())), len(set(non_repeatregion_vcf_dict.keys()) & set(non_repeatregion_simu_dict.keys()))
	pickle_file = '/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record_simu_picklefile.txt'%(length, depth)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = ( simu_dict, repeatregion_simu_dict, non_repeatregion_simu_dict)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
def main():
	task_pool = Pool(processes = 8)
	for length in [75, 150]:
		for depth in [10, 30, 50]:
			#'''
			print "processing length:%s, depth:%s"%(length, depth)
			pjobs_ids = task_pool.apply_async(get_simu_dict, args=(length, depth,))
			#get_simu_dict(length, depth)
	print "Waiting for all subprocesses done..."
	task_pool.close()
	task_pool.join()
	
	print 'All subprocesses done.'
	
	
main()
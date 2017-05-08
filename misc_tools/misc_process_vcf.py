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
def get_vcf_dict(length, depth):
	print "processing length:%s, depth:%s"%(length, depth)
	repeat_region = csv.reader(open("/zzh_gpfs02/zhangyanfang/repeat_region_analysis/all_repeat_region_interval.txt", "rU"), delimiter = "\t")
	repeat_region_dict = {}
	for index, line in enumerate(repeat_region):
		repeat_region_dict.setdefault(line[0], []).append((line[1], line[2]))


	vcf_file = "/zzh_gpfs02/yanmingchen/wholegenome/hapler_combine_order/simu_%s_%sX.filter.uniq.concordant.merge.sort.dedup.AddGroup.snp.haplo.g.vcf.haplercaller.vcf"%(length, depth)
	handle = csv.reader(open(vcf_file, "rU"), delimiter = "\t")
	vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict, haplo_type_set = {},{},{}, []
	haplo_type_r = []
	haplodict = {"1/1": "-",  "0/1": "+", "1/0": "+", "0/0": "-"}
	for index, line in enumerate(handle):
		if index % 10000 == 0:
			print "Processed 10000 lines in simu %s %s."%(length, depth)
		if "chr" in line[0] and len(line) == 10:
			chromsome = line[0]
			position  = int(line[1])
			ref_nucle = line[3]
			snv_nucle = line[4]
			haplo_type_r.append(line[9].split(":")[0])
			in_repeat_region = 0
			try:
				for region in repeat_region_dict[chromsome]:
					if int(region[0]) <= position-1 <= int(region[1]):
						in_repeat_region = 1
			except KeyError:
				pass
			try:
				haplo_type = haplodict[line[9].split(":")[0]]
			except KeyError:
				haplo_type = "-"
			#print chromsome, ref_nucle, snv_nucle, haplo_type
			vcf_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
			if in_repeat_region == 1:
				repeatregion_vcf_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
			else:
				non_repeatregion_vcf_dict[(chromsome, position)] = (ref_nucle, snv_nucle, haplo_type, in_repeat_region)
			haplo_type_set.append(haplo_type)
			#sys.exit(0)
	#'''
	print Counter(haplo_type_set), Counter(haplo_type_r)




	#How many position be detected
	#print len(set(vcf_dict.keys())), len(set(simu_dict.keys())), len(set(vcf_dict.keys()) & set(simu_dict.keys()))
	#print len(set(repeatregion_vcf_dict.keys())), len(set(repeatregion_simu_dict.keys())), len(set(repeatregion_vcf_dict.keys()) & set(repeatregion_simu_dict.keys()))
	#print len(set(non_repeatregion_vcf_dict.keys())), len(set(non_repeatregion_simu_dict.keys())), len(set(non_repeatregion_vcf_dict.keys()) & set(non_repeatregion_simu_dict.keys()))
	pickle_file = '/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record_vcf_picklefile.txt'%(length, depth)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = ( vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
def clean_dict(dict1):
	for (key, value) in dict1.items():
		if len(key.split("_")) == 1: 
			pass
def caculate_dict( simu_dict, vcf_dict):
	Simu_number = len(simu_dict)
	vcf_number  = len(vcf_dict)
	com_pos     = set(vcf_dict.keys()) & set(simu_dict.keys())
	PP          = len(set(vcf_dict.keys()) & set(simu_dict.keys()))
	PPR         = float(PP) / float(vcf_number)
	N_position  = 0
	for pos in com_pos:
		#print vcf_dict[pos], simu_dict[pos]
		
		#if len(simu_dict[pos][1]) == 2:
		simurecord = (simu_dict[pos][0], [x for x in simu_dict[pos][1] if x != simu_dict[pos][0]][0], simu_dict[pos][2], simu_dict[pos][3])
		#print vcf_dict[pos], simurecord, simu_dict[pos]
		#sys.exit(0)
		if vcf_dict[pos] == simurecord:
			N_position += 1
	N_Accuracy	= N_position
	Accuracy_Rate = float(N_Accuracy) / float(Simu_number)
	Positive_rate  = float(N_position) / float(vcf_number)
	N_FP = PP - N_position + len(set(vcf_dict.keys())) - PP
	FPR = float(N_FP) / float(vcf_number)
	N_FN = Simu_number - PP
	FNR  = float(Simu_number - PP )/ float(Simu_number)
	return  [ Simu_number, vcf_number, N_Accuracy, Accuracy_Rate, PP, PPR, N_position, Positive_rate, N_FP, FPR, N_FN, FNR]
def caculate_info(index, length, depth):
	vcf_file = '/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record_vcf_picklefile.txt'%(length, depth)
	simu_file = '/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record_simu_picklefile.txt'%(length, depth)
	output_file = '/zzh_gpfs02/yanmingchen/WGS_simulation/test/simu_%s_%sX_record_statics.txt'%(length, depth)
	pickle_tuple = pickle.load(open(vcf_file, 'rb'))
	vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict = pickle_tuple[0], pickle_tuple[1], pickle_tuple[2]
	pickle_tuple = pickle.load(open(simu_file, 'rb'))
	simu_dict, repeatregion_simu_dict, non_repeatregion_simu_dict = pickle_tuple[0], pickle_tuple[1], pickle_tuple[2]
	output_handle = csv.writer(open(output_file, "w"), delimiter = "\t")
	output_handle.writerow(["Sample", "Length", "Depth", "#Simu", "#Find", "Accuracy", "Accuracy_Rate", "#PP", "PPR", "#Positive", "Positive_Rate", "#FP","FPR", "#FN", "FNR"])
	output_handle.writerow([index, length, depth ] + caculate_dict( simu_dict, vcf_dict))
	output_handle.writerow([index, length, depth ] + caculate_dict( repeatregion_simu_dict, repeatregion_vcf_dict))
	output_handle.writerow([index, length, depth ] + caculate_dict(non_repeatregion_simu_dict, non_repeatregion_vcf_dict))
def main():
	task_pool = Pool(processes = 8)
	index = 0
	for length in [75, 150]:
		for depth in [10, 30, 50]:
			#'''
			index += 1
			#pjobs_ids = task_pool.apply_async(get_vcf_dict, args=( length, depth,))
			try:
				caculate_info(index, length, depth)
			except IOError:
				pass
	print "Waiting for all subprocesses done..."
	task_pool.close()
	task_pool.join()
	
	print 'All subprocesses done.'
	
			
main()
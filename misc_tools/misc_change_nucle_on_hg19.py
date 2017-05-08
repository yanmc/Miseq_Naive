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
def exchange_nucl(key, value, gennome_record):
	
	my_seq = gennome_record.seq.tomutable()
	for info in value:
		my_seq[info[0]] = info[1]
	return SeqRecord(my_seq.toseq(), id = key, description = '')
def main():
	snv_records_file = "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/one_snv_analysis/new/all_snv_from_extend.vcf"
	genome_file = "/zzh_gpfs/data/ref_genome/hg19_genome/hg19.fa"
	gennome_record_handle = SeqIO.parse(open(genome_file, "rU"), "fasta")
	gennome_record = {}
	for record in gennome_record_handle:
		gennome_record[record.id] = record
	handle = csv.reader(open(snv_records_file, "rU"), delimiter = "\t")
	handle.next()
	handle.next()
	info, result = [], {}
	for index, line in enumerate(handle):
		if index %1000 == 0:
			print "%s has Done!"%index
		query_chr = line[0]
		query_position = int(line[4]) -1
		subject_chr = line[6]
		subject_position =  int(line[9]) -1
		ref_nucle = line[10]
		alt_nucle = line[11]
		#info = [query_chr, query_position, subject_chr, subject_position, ref_nucle, alt_nucle]
		info.append((query_chr, query_position, alt_nucle))
		info.append((subject_chr, subject_position, ref_nucle))
		#print info
		#print  [x for x in gennome_record.keys()]
		#print gennome_record[info[0]].seq[info[1]], gennome_record[info[2]].seq[info[3]]
		#print  gennome_record[info[0]].seq[info[1] -5: info[1] + 6], gennome_record[info[2]].seq[info[3] - 5 : info[3] +6]
		#gennome_record = exchange_nucl(info, gennome_record)
		#print gennome_record[info[0]].seq[info[1]], gennome_record[info[2]].seq[info[3]]
		#if index == 10:
			#sys.exit(0)
			#break
	#print info
	for i, change_record in groupby(info, lambda s:s[0]):
		result.setdefault(i, []).append([v[1:] for v in change_record][0])
	#print result
	os.system("rm ./changed_hg19.fasta")
	outfile = open("./changed_hg19.fasta", "wa+")
	for key, value in result.items():
		print "Change %s ..."%key
		my_gennome_record = exchange_nucl(key, value, gennome_record[key])
		SeqIO.write(my_gennome_record, outfile, "fasta")
		
main()
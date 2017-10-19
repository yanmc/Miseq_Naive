#!/usr/bin/env python
# encoding: utf-8
"""
2.0-process-igblast-output.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
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
try:
    import cPickle as pickle
except ImportError:
    import pickle
import math
from mpl_toolkits.mplot3d import Axes3D
import copy

def get_germline_gene(germline_type, chain_type):
	if germline_type == 'V' and chain_type == "H":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGHV']
	elif germline_type == 'J' and chain_type == "H":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGHJ']
	
	elif germline_type == 'V' and chain_type == "K":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGKV']
	elif germline_type == 'J' and chain_type == "K":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGKJ']
	
	elif germline_type == 'V' and chain_type == "L":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGLV']
	elif germline_type == 'J' and chain_type == "L":
		germline_gene_list = HUMAN_GERMLINE['HUMANIGLJ']
	return germline_gene_list

def get_all_gene_name(germline_gene_list):
	return sorted(set([x.split('*')[0] for x in germline_gene_list if "OR" not in x]))
	
def get_VJ_recombanation(line):
	Sample_ID = line[0][1:]
	VJ_pair   = tuple([line[1].split(',')[0].split('*')[0], line[3].split(',')[0].split('*')[0]])
	return Sample_ID, VJ_pair
def main():
	print "Begin!"
	V_genes = get_all_gene_name(get_germline_gene("V", chain_type))
	J_genes = get_all_gene_name(get_germline_gene("J", chain_type))
	print V_genes, J_genes
	real_reads_infile = "%s/%s_H_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name)
	real_reads_list = []
	for record in SeqIO.parse(real_reads_infile, "fasta"):
		real_reads_list.append(record.id)
	print "real_reads_list", len(real_reads_list), real_reads_list[0]
	infile = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	print "Processing %s " %infile
	handle = open(infile,"rU")
	reader = csv.reader(handle,delimiter = "\t")
	Recombanation_dict = {}
	for line in reader:
		#print line
		if len(line) > 1 and line[-2] == "Yes":
			Sample_ID, VJ_pair = get_VJ_recombanation(line)
			#print list(Sample_ID), VJ_pair
			if Sample_ID in real_reads_list:
				Recombanation_dict[Sample_ID] = VJ_pair
	#Processing DNA
	real_reads_CDR3_infile = "%s/%s_H_real_reads_CDR3.fasta"%(prj_tree.reads, prj_name)
	#real_reads_CDR3_dict = SeqIO.index(real_reads_CDR3_infile, "fasta")
	real_reads_CDR3_dict = {}
	for record in SeqIO.parse(real_reads_CDR3_infile, "fasta"):
		#if str(record.seq) != "GCGAGAAGGGCCTTAGCAGCAGCTGGTATCCCGGTCAGGGGCTGGTTCGACCCCTGGGGCCAGGGA" and str(record.seq) != "GCAAGAGCAGCAGAAGGGCGGGGGGGCTACTACTACTACTACATGGACGTCTGGGGCAAAGGG":
		#	real_reads_CDR3_dict[record.id] = record
		real_reads_CDR3_dict[record.id] = record
	real_reads_CDR3_record = "%s/%s_H_real_reads_CDR3_unique.record"%(prj_tree.reads, prj_name)
	clone_dict = {}
	for index, line in enumerate(csv.reader(open(real_reads_CDR3_record, "rU"), delimiter = "\t")):
		for sample_ID in line:
			try:
				CDR3_seq = real_reads_CDR3_dict[sample_ID].seq#.translate()
				VJ_pair = Recombanation_dict[sample_ID]
				clone_dict.setdefault(VJ_pair + (CDR3_seq,), []).append(sample_ID)
			except KeyError:
				pass
	
	writer = csv.writer(open("%s/%s_H_clone_frequency_all2.txt"%(prj_tree.data, prj_name), "w"), delimiter = "\t")
	writer.writerow(["V gene", "J gene", "Nucl Seq", "AA Seq", "Frequency (%)", "Number"])
	for index, (key, value) in enumerate(sorted(clone_dict.items(), key=lambda d: len(d[1]), reverse = True)):
		result =  [key[0], key[1], str(key[2]), str(key[2].translate()), (float(len(value))/len(real_reads_CDR3_dict)) *100 , len(value)]
		
		#if index < 20:
		writer.writerow(result)
			#sys.exit(0)
	
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	print time.time()
	start = time.time()
	'''
	if "K" in prj_name:
		chain_type = "K"
		main()
	elif "L" in prj_name:
		chain_type = "L"
		main()
	elif "H" in prj_name:
		chain_type = "H"
		main()
	else:
		#for chain_type in ("H", "K", "L"):
		for chain_type in ("H"):
			print "Processing %s chain..."%chain_type
			main()
	'''
	for chain_type in ("H"):
		print "Processing %s chain..."%chain_type
		main()
	end = time.time()
	print prj_name, end-start
	print "Finished"
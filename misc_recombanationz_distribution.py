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
				Recombanation_dict.setdefault(VJ_pair, []).append(Sample_ID)
	#print "Recombanation_dict", Recombanation_dict
	#for index, (key, value) in enumerate(sorted(Recombanation_dict.items())):
	#	print key, value
	#	sys.exit(0)
	#	if key[0] in V_genes and  key[1] in J_genes:
	#		print key, value
	#X, Y = V_genes, J_genes
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	'''
	y = []
	for J_gene in J_genes:
		y.append(tuple(copy.deepcopy(range(1, len(V_genes)+1))))
	x = range(1, len(J_genes)+1)
	print y , x
	for xe, ye in zip(x, y):
		   plt.scatter([xe] * len(ye), ye)
	
	plt.xticks(range(1, len(V_genes)+1))
	plt.axes().set_xticklabels(V_genes)
	plt.yticks(range(1, len(J_genes)+1))
	plt.axes().set_yticklabels(J_genes)
	'''
	#'''
	#X_axis, Y_axis = list(np.array(range(1,len(V_genes)+1))), list(np.array(range(1,len(V_genes)+1)))
	plt.figure(figsize=(33, 10))
	colors = "red"
	area = np.zeros( (len(J_genes), len(V_genes)) )
	
	for J_index, J_gene in enumerate(J_genes):
		for V_index, V_gene in enumerate(V_genes):
			try:
				#pair_num = len(Recombanation_dict[tuple([V_gene, J_gene])])
				pair_num = float(len(Recombanation_dict[tuple([V_gene, J_gene])]))/len(real_reads_list) 
			except KeyError:
				pair_num = 0
				#pair_num = 0
			#area[J_index, V_index] = np.pi * (1 * math.log(pair_num, 2))**2
			#area[J_index, V_index] =  math.log(pair_num, 2)
			#area[J_index, V_index] = np.pi * (pair_num)**2
			area[J_index, V_index] = pair_num *600
	for index, line in enumerate(area):
		
		Y_axis = [index for x in range(1,len(V_genes)+1)]
		X_axis = range(1,len(V_genes)+1)
		#print help()
		plt.scatter(X_axis, Y_axis, s=line, c=colors, alpha=0.5)
	
	plt.xticks(range(1, len(V_genes)+1))
	plt.axes().set_xticklabels(V_genes, rotation = 45, ha = 'right')#, multialignment = 'right')
	plt.yticks(range(0, len(J_genes)))
	plt.axes().set_yticklabels(J_genes)
	plt.axes().set_xlim(0,86)
	plt.grid(True)
	#plt.axes().set_ylim(0,7)
	
	#plt.xticks( [x/len(V_genes) for x in X_axis], V_genes)
	#plt.yticks( [x/len(J_genes) for x in Y_axis], J_genes)
	#'''
	
	'''
	print X_axis, Y_axis, area
	c_dict = {1:"r", 2:"g", 3:"b", 4:"c", 5:"y", 6:"k"}
	for index, line in enumerate(area):
		xs = np.arange(len(line))
		cs = c_dict[index+1] * len(xs)
		ax.bar(xs, line, zs=index, zdir='y', color=cs, alpha=0.8)
	
	x_ticks, x_ticks_label = xs, V_genes
	ax.set_xticklabels(x_ticks_label)
	#ax.set_yticks(x_ticks, x_ticks_label)
	#ax.set_zticks(x_ticks, x_ticks_label)
	ax.tick_params
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	os.system("rm %s/%s_get_recombanation_info.png"%(prj_tree.figure, prj_name))
	'''
	plt.savefig("%s/%s_get_recombanation_info.png"%(prj_tree.figure, prj_name), dpi = 300)

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
	
	#infiles = glob.glob("%s/*_real_reads_Variable_region*.fasta"%prj_tree.reads)
	#for infile in infiles:
	#	tanslate_file(infile)
	
	#infiles = glob.glob("%s/*_real_reads_CDR3*.fasta"%prj_tree.reads)
	#for infile in infiles:
	#	tanslate_file(infile)
	
	
	
	end = time.time()
	print prj_name, end-start
	print "Finished"
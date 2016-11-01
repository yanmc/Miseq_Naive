#!/usr/bin/env python
# encoding: utf-8
"""
misc-merge_IgD_IgM_data.py 

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

def get_all_gene_name(germline_gene_list):
	return sorted(set([x.split('*')[0] for x in germline_gene_list]))

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
def get_naive_right_allele(data_dict, pic_type, germline_gene_list, germline_type):
	gene_number_dict, result_dict = {}, {}
	total_number = sum(data_dict.values())
	for (key, value) in data_dict.items():
		gene_number_dict.setdefault(key.split('*')[0], []).append((key, value))
	for (key, value) in  gene_number_dict.items():
		print key, value
		if len(value) > 1:
			sorted_value = sorted(value, key = lambda z : z[1], reverse=True)
			if float(sorted_value[1][1])/float(sorted_value[0][1]) < 0.001:
				sorted_value = [sorted_value[0][0]]
			else:
				sorted_value = [x[0] for x in sorted_value[:2]]
		else:
			sorted_value = [value[0][0]]
		result_dict[key] = sorted_value
	print result_dict
	result_file = csv.writer(open('/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s_%s_%s_gene_usage%s.txt'%(prj_name, chain_type, germline_type, pic_type),"w"), delimiter="\t")
	result_file.writerow(["Germline gene", "Allele 1", "Allele 2"])
	for (key, value) in result_dict.items():
		print key, value
		if len(value) == 1:
			result_file.writerow([key, value[0].split('*')[1], 0])
		else:
			result_file.writerow([key, value[0].split('*')[1], value[1].split('*')[1]])
def plot_gene_usage(result_dict, pic_type, germline_gene_list, germline_type):
	SBG = StackedBarGrapher()
	gene_number_dict = {}
	total_number = sum(result_dict.values())
	for (key, value) in result_dict.items():
		gene_number_dict.setdefault(key.split('*')[0], []).append(float(value)/total_number*100)
	max_len_value = 0
	for key, value in  gene_number_dict.items():
		if max_len_value < len(value):
			max_len_value = len(value)
	gene_names = get_all_gene_name(germline_gene_list)
	data = np.zeros((len(gene_names),max_len_value))
	for index, gene in enumerate(gene_names):
		try:
			gene_number = gene_number_dict[gene]
			print gene_number
			for n_index, number in enumerate(sorted(gene_number, reverse=True)):
				
				data[index][n_index] = number
		except KeyError:
			pass		

	d_labels = gene_names
	d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777', '#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
	rows = ['%d%%' % x for x in range(75, 101)]
	#rows.insert(0, "<75%")
	colors = plt.cm.rainbow(np.linspace(0, 1, 10))
	#d_colors = colors
	fig = plt.figure()
	max_y_tick = math.ceil(max(np.sum(data, axis=1)))
	y_ticks_tick, y_ticks_label, tick_log = [0],[0], 0
	
	while max_y_tick % 2 != 0:
		max_y_tick += 1
	while tick_log <= max_y_tick:
		tick_log += 2
		y_ticks_tick.append(tick_log*1)
		y_ticks_label.append(tick_log*1)
	while len(y_ticks_tick) > 14:
		if len(y_ticks_tick) %2 == 0:
			y_ticks_tick, y_ticks_label = y_ticks_tick[1:][::2], y_ticks_label[1:][::2]
			y_ticks_tick.insert(0,0), y_ticks_label.insert(0,0)
		else:
			y_ticks_tick, y_ticks_label = y_ticks_tick[::2], y_ticks_label[::2]
	y_ticks = [y_ticks_tick, y_ticks_label]
	ax = fig.add_subplot(111)
	SBG.stackedBarPlot(ax,
	                   data,
	                   d_colors,
	                   xLabels=d_labels,
                       edgeCols=['#000000']*len(d_colors),
	                   yTicks= y_ticks,
					   ylabel = 'Percent of reads (%)',
					   gap =.2,
                       endGaps=True
	                  )
	
	plt.title("%s %s %s gene usage"%(prj_name, chain_type, germline_type))

	fig.subplots_adjust(bottom=0.4)
	plt.tight_layout()
	for t in ax.xaxis.get_ticklabels():
		t.set_horizontalalignment('center')
		if str(t).split("'")[1] in list(gene_number_dict.keys()):
			t.set_color('blue')
			
	plt.savefig('/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s_%s_%s_gene_usage%s.png'%(prj_name, chain_type, germline_type, pic_type), dpi=300)
	del fig
	plt.close()
def main():
	germline_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta")
	
	for germline_type in ('V','J'):
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		pickle_file_IgM = '/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s_IgM/tmp/%s_IgM_gene_usage_info_dump_%s_H%s'%(prj_name, prj_name, germline_type, pic_type)
		pickle_file_IgD = '/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s_IgD/tmp/%s_IgD_gene_usage_info_dump_%s_H%s'%(prj_name, prj_name, germline_type, pic_type)
		f_IgM = open(pickle_file_IgM, 'rb')
		pickle_tuple_IgM = pickle.load(f_IgM)
		f_IgD = open(pickle_file_IgD, 'rb')
		pickle_tuple_IgD = pickle.load(f_IgD)
		columns_IgM, geneusage_data_IgM, gene_usage_ids_IgM = pickle_tuple_IgM[0], pickle_tuple_IgM[1], pickle_tuple_IgM[2]
		columns_IgD, geneusage_data_IgD, gene_usage_ids_IgD = pickle_tuple_IgD[0], pickle_tuple_IgD[1], pickle_tuple_IgD[2]
		IgM_geneusage_dict, IgD_geneusage_dict = {}, {}
		for i in range(len( pickle_tuple_IgM[0])):
			IgM_geneusage_dict[columns_IgM[i]] = sum(geneusage_data_IgM[i])
		for i in range(len( pickle_tuple_IgD[0])):
			IgD_geneusage_dict[columns_IgD[i]] = sum(geneusage_data_IgD[i])
		#germline_ids = set(germline_fasta.keys())
		germline_ids = set(IgM_geneusage_dict) | set(IgD_geneusage_dict)
		result_dict = {}
		for germline_id in germline_ids:
			try:
				IgM_number = IgM_geneusage_dict[germline_id]
			except KeyError:
				IgM_number = 0
			try:
				IgD_number = IgD_geneusage_dict[germline_id]
			except KeyError:
				IgD_number = 0
			print germline_id, IgM_number, IgD_number
			result_dict[germline_id] =  IgM_number + IgD_number
		plot_gene_usage(result_dict, pic_type, germline_gene_list, germline_type)
		get_naive_right_allele(result_dict, pic_type, germline_gene_list, germline_type)
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	pic_type = "_unique"
	chain_type = "H"
	start = time.time()
	for prj_name in ["Naive", "Nsw"]:
		main()
	
	end = time.time()
	print prj_name, end-start
	print "Finished"
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
	result_file = csv.writer(open('/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s/%s_%s_%s_gene_usage%s.txt'%(prj_name, prj_name, chain_type, germline_type, pic_type),"w"), delimiter="\t")
	result_file.writerow(["Germline gene", "Allele 1", "Allele 2"])
	max_freq_allele_dict = {}
	for (key, value) in result_dict.items():
		print key, value
		if len(value) == 1:
			result_file.writerow([key, value[0].split('*')[1], 0])
			max_freq_allele_dict[key] = value[0]
		else:
			result_file.writerow([key, value[0].split('*')[1], value[1].split('*')[1]])
			max_freq_allele_dict[key] = value[0]
	return max_freq_allele_dict
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
			
	plt.savefig('/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s/%s_%s_%s_gene_usage%s.png'%(prj_name, prj_name, chain_type, germline_type, pic_type), dpi=300)
	del fig
	plt.close()
def get_novel_allele(result_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta):
	gene_reads_id_dict = {}
	for (key, value) in result_ids_dict.items():
		gene_reads_id_dict.setdefault(key.split('*')[0], []).extend(value)
	for key,value in gene_reads_id_dict.items():
		print key, type(value), len(value)
		sys.exit(0)
	for genmline_gene in germline_gene_list:
		try:
			max_freq_allele = max_freq_allele_dict[genmline_gene]
		except KeyError:
			continue
		all_reads_ids_list = gene_reads_id_dict[genmline_gene]
		get_mutation_patterns(all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta)
		get_mutation_spectrum(all_reads_ids_list, max_freq_allele)

def get_mutation_patterns(all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta):
	mutation_patterns_dict = {}
	print type(all_reads_ids_list), len(all_reads_ids_list)
	for read_id in all_reads_ids_list:
		ref_seq_id, test_seq_id = max_freq_allele, read_id
		print ref_seq_id, type(test_seq_id),len(test_seq_id), type(germline_fasta), type(unique_real_reads_fasta)
		ref_seqrecord, test_seqrecord = germline_fasta[ref_seq_id], unique_real_reads_fasta[test_seq_id]
		ref_seq_id, test_seq_id = ref_seq_id.replace('/','').replace('*',''), test_seq_id.replace('/','').replace('*','')
		out = open('%s_%s_pair.fasta'%(test_seq_id,ref_seq_id),'w')
		SeqIO.write(ref_seqrecord, out, 'fasta')
		SeqIO.write(test_seqrecord ,out, 'fasta')
		out.close()
		#my_ref_len = len(rank_germ[i])
		file_for_clustalw = '%s_%s_pair.fasta'%(test_seq_id,ref_seq_id)
		do_clustalw(file_for_clustalw)
		clustalw_result = '%s_%s_pair.aln'%(test_seq_id,ref_seq_id)
		mutation_patterns_dict = caculate_mutation_patterns(clustalw_result, read_id, mutation_patterns_dict)
		os.system("rm %s_%s_pair.fasta"%(test_seq_id,ref_seq_id))
		os.system("rm %s_%s_pair.aln"%(test_seq_id,ref_seq_id))
		os.system("rm %s_%s_pair.dnd"%(test_seq_id,ref_seq_id))
	mutation_patterns_group = {}
	for (key, value) in mutation_patterns_dict:
		mutation_patterns_group.setdefault(value[0], []).append((key, value[0]))
	data = np.zeros(len(mutation_patterns_group), len(germline_fasta[max_freq_allele]))
	for index, (group_number, value) in enumerate(mutation_patterns_group):
		for item in value:
			position = item[1][0]
			data[index][position] += 1
	ref_seq_id_name = ref_seq_id.split('*')[0]
	mutation_patterns_file = open('/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/%s/%s_%s_mutation_patterns.txt'%(prj_name, prj_name, ref_seq_id_name), 'w')
	print data
	mutation_patterns_writer = csv.writerows(data)
	mutation_patterns_file.close()
def get_mutation_spectrum(all_reads_ids_list, max_freq_allele, germline_fasta):
	pass

def remove_ref_insertion(ref, tst):
	while ref.find("-") >= 0:
		index = ref.index("-")
		ref = ref[ :index] + ref[index + 1 :]
		tst = tst[ :index] + tst[index + 1 :]
	trimmed_ref_len = len(ref)
	return ref, tst, trimmed_ref_len
def parse_pair_clustal2(align_file):
	""" parse paired alignment file """
	alignment 	= AlignIO.read(align_file, "clustal")
	seq1, seq2 	= alignment[0].seq.tostring(), alignment[1].seq.tostring()
	seq1, seq2, trimmed_ref_len = remove_ref_insertion(seq1, seq2)
	zip_seqs 	= zip(seq1, seq2)
	#zip_seqs	= [(x, y) for x, y in zip_seqs if x !="-" and y != "-"]		# remove insertions from either one
	mismatches	= sum([x != y for x, y in zip_seqs])
	mutation_patterns_result = []						# count total matches
	for index, (x, y) in enumerate(zip_seqs):
		if x != y:
			position = index + 1
			mutation_patterns_result.append((position, x, y))
	return mutation_patterns_result, mismatches
def caculate_mutation_patterns(clustalw_result, read_id, mutation_patterns_dict):
	result = csv.writer(open('%s_mutation_patterns.txt'%prj_name,'a+'),delimiter = '\t')
	fs = glob.glob(clustalw_result)
	for infile in fs:
		print "processing %s"%infile
		mutation_patterns_result, mismatches = parse_pair_clustal2(infile)
		mutation_patterns_dict[read_id] = (mismatches, mutation_patterns_result)
	return mutation_patterns_dict
def do_clustalw(file_for_clustalw):
	infiles = glob.glob(file_for_clustalw)

	#clustalw_exe = r"/Applications/clustalw-2.1-macosx/clustalw2"
	clustalw_exe = r"/zzh_gpfs/apps/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
	assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
	for in_file in infiles:
		print "Processing %s ......."%in_file
		in_file = in_file.replace('&','\&')
		in_file = in_file.replace('*','\*')
		clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file)
		stdout, stderr = clustalw_cline()

def main():
	os.system("mkdir ../%s"%prj_name)
	germline_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta")
	IgD_unique_real_reads_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgD/reads/Naive_IgD_H_real_reads_V_region_unique.fasta")
	IgM_unique_real_reads_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/reads/Naive_IgM_H_real_reads_V_region_unique.fasta")
	unique_real_reads_fasta =dict(IgD_unique_real_reads_fasta, **IgM_unique_real_reads_fasta)
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
		IgM_geneusage_dict, IgD_geneusage_dict, IgM_geneusage_ids_dict, IgD_geneusage_ids_dict = {}, {}, {}, {}
		for i in range(len( pickle_tuple_IgM[0])):
			IgM_geneusage_dict[columns_IgM[i]] = sum(geneusage_data_IgM[i])
			IgM_geneusage_ids_dict[columns_IgM[i]] = (gene_usage_ids_IgM[i][0]) | gene_usage_ids_IgM[i][1]
		for i in range(len( pickle_tuple_IgD[0])):
			IgD_geneusage_dict[columns_IgD[i]] = sum(geneusage_data_IgD[i])
			IgD_geneusage_ids_dict[columns_IgD[i]] = gene_usage_ids_IgD[i][0] | gene_usage_ids_IgD[i][1]
		#germline_ids = set(germline_fasta.keys())
		germline_ids = set(IgM_geneusage_dict) | set(IgD_geneusage_dict)
		result_dict,  result_ids_dict = {}, {}
		for germline_id in germline_ids:
			try:
				IgM_number = IgM_geneusage_dict[germline_id]
				IgM_reads_id = IgM_geneusage_ids_dict[germline_id]
			except KeyError:
				IgM_number = 0
				IgM_reads_id = []
			try:
				IgD_number = IgD_geneusage_dict[germline_id]
				IgD_reads_id = IgD_geneusage_ids_dict[germline_id]
			except KeyError:
				IgD_number = 0
				IgD_reads_id = []
			result_dict[germline_id] =  IgM_number + IgD_number
			result_ids_dict[germline_id] = IgM_reads_id + IgD_reads_id
		plot_gene_usage(result_dict, pic_type, germline_gene_list, germline_type)
		max_freq_allele_dict = get_naive_right_allele(result_dict, pic_type, germline_gene_list, germline_type)
		print "get_novel_allele"
		for key,value in result_ids_dict.items():
			print key, type(value), len(value)
			sys.exit(0)
		get_novel_allele(result_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta)
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
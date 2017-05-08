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

def process_recomb(infile, germline_type, germline_gene_list):
	real_reads_list = []
	print "Processing %s " %infile
	csv_handle = open(infile,"rU")
	reader = csv.reader(csv_handle,delimiter = "\t")
	total, count = 0, 0
	for index, line in enumerate(reader):
		if len(line) != 1:
			count += 1
		recomb_result = Recomb(line)
		if recomb_result.v == "N/A" or recomb_result.j == "N/A":
			recomb_result.set_cover_vj("No")
		else:
			recomb_result.set_cover_vj("Yes")
		#print recomb_result.v, recomb_result.d, recomb_result.j, recomb_result.productive, recomb_result.cover_vj
		if germline_type == "V":
			if recomb_result.productive == "Yes" and recomb_result.cover_vj == "Yes" and recomb_result.v.split('*')[0] in germline_gene_list:
				real_reads_list.append(recomb_result.qid)
		if germline_type == "J":
			if recomb_result.productive == "Yes" and recomb_result.cover_vj == "Yes" and recomb_result.j.split('*')[0] in germline_gene_list:
				real_reads_list.append(recomb_result.qid)
	csv_handle.close()
	return real_reads_list, count, index+1

def process_alignment(infile, germline_type, real_reads_list, coverage_cutoff):
	#NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = {}, {}, [], []
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = {}, {}, {}
	print "Processing %s " %infile
	handle = open(infile,"rU")
	reader = csv.reader(handle,delimiter = "\t")
	count1,count2,count3 = 0, 0, 0
	for line in reader:
		assign_result = MyAlignment(line) #Use my_tools's MyAlignment CLASS
		coverage_rate = assign_result.coverage_rate
		if assign_result.assign_type == germline_type:
			#if total % 1000 == 0:
			#	print "%s Done!"%total
			if assign_result.identity == float(100.00) and assign_result.sstart == 1 and coverage_rate >= coverage_cutoff: #and assign_result.send == assign_result.slen:
				count1 +=1
				#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
				#NoMut_reads_list_all.append(assign_result.qid)
				#NoMut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#NoMut_reads_list.append(assign_result.qid)
					NoMut_alignment_dict.setdefault(assign_result.sid, []).append(assign_result.qid)
					maturation_rate_dict.setdefault(assign_result.identity, []).append(assign_result.qid)
			elif coverage_rate >= coverage_cutoff:
				count2 +=1
				#Mut_reads_list_all.append(assign_result.qid)
				#Mut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#Mut_reads_list.append(assign_result.qid)
					Mut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
					maturation_rate_dict.setdefault(assign_result.identity, []).append(assign_result.qid)
			else:
				count3 +=1
		#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)	
	#print total,len(NoMut_reads_list), len(Mut_reads_list), count
	handle.close()
	return NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict

def def_score_function(alignment_dict):
	subject_id_value_dict, value_list = {}, []
	for value in alignment_dict.values():
		#print value
		for item in value:
			#processed_item = (item[0],float(1)/len(value))
			#processed_item = (item[0],1)#simple way
			processed_item = (item[0],float(item[1])/len(value))
			#processed_item = (item[0],float(item[1]))
			value_list.append(processed_item)
	return value_list


def calculate_score(value_list):
	unique_ID = []
	for item in value_list:
		flag=1
		for i in range(0,len(unique_ID)):
			if item[0] == unique_ID[i][0]:
				unique_ID[i][1] += item[1]
				flag=0
				break
		if flag==1:
			unique_ID.append(list(item))
	return unique_ID
def get_germ_ids(x, y):
	return set(x)|set(y)

def get_alignment_info(IGBLAST_assignment_file, germline_type, real_reads_list, coverage_cutoff):
	index = IGBLAST_assignment_file.split('/')[-1].split('_')[2]
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list, coverage_cutoff)
	pickle_file = '%s/%s_get_assignment_info_dump_%s_%s_%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, index)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = (NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()

def process_dump(infile):
	f = open(infile, 'rb')
	pickle_tuple = pickle.load(f)
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = pickle_tuple[0], pickle_tuple[1], pickle_tuple[2]
	f.close()
	return	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict

def combine_dict_and_list(NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict, NoMut_alignment_dict_part, Mut_alignment_dict_part, maturation_rate_dict_part):
	
	NoMut_alignment_dict_combined = combine_dict_part(NoMut_alignment_dict, NoMut_alignment_dict_part)
	Mut_alignment_dict_combined = combine_dict_part(Mut_alignment_dict, Mut_alignment_dict_part)
	maturation_rate_dict_combined = combine_dict_part(maturation_rate_dict, maturation_rate_dict_part)
	
	return NoMut_alignment_dict_combined, Mut_alignment_dict_combined, maturation_rate_dict_combined

def combine_dict_part(Mut_alignment_dict, Mut_alignment_dict_part):
	Mut_alignment_dict_combined = {}
	germline_ids = get_germ_ids(Mut_alignment_dict.keys(), Mut_alignment_dict_part.keys())
	for germline_id in germline_ids:
		try:
			Mut_reads = Mut_alignment_dict[germline_id]
		except KeyError:
			Mut_reads = []
		try:
			Mut_reads_part = Mut_alignment_dict_part[germline_id]
		except KeyError:
			Mut_reads_part = []
		Mut_alignment_dict_combined[germline_id] =  Mut_reads + Mut_reads_part
	return Mut_alignment_dict_combined
def continuous_subtraction_list(Mut_reads_num_i):
	num = 0
	for item in Mut_reads_num_i:
		item[1] = item[1] - num
		num = num + item[1]
	return [x[1] for x in Mut_reads_num_i]

def get_identity_list(Mut_reads, identity_dict):
	Mut_reads_num_i = [[x, 0] for x in range(75, 101)]
	for read_id in Mut_reads:
		identity = identity_dict[read_id]
		for i in range(75,101):	
			if identity < i:
				Mut_reads_num_i[i - 75][1] += 1
	return Mut_reads_num_i

def get_max_y_tick(max_value):
	#return int(math.ceil(float(max_value)/10000))
	max_part = int(('').join([x for x in str(max_value) if x != "0"]))
	max_part_log = int(max_value/max_part)
	while max_part % 2 != 0:
		max_part += 1
	return max_part, max_part_log
def clear_VJcoverage80_reads(alignment_dict, real_reads):
	for key, value in alignment_dict.items():
		new_value = set(real_reads) & set(value)
		if len(new_value) == 0:
			del alignment_dict[key]
		else:
			alignment_dict[key] = new_value
	return alignment_dict

def get_dump_data(germline_type):
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict= {}, {}, {}
	get_assignment_info_dump_files = glob.glob('%s/%s_get_assignment_info_dump_%s_%s_*'%(prj_tree.tmp, prj_name, germline_type, chain_type))
	for infile in get_assignment_info_dump_files:
		NoMut_alignment_dict_part, Mut_alignment_dict_part, maturation_rate_dict_part = process_dump(infile)
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = combine_dict_and_list(NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict, NoMut_alignment_dict_part, Mut_alignment_dict_part, maturation_rate_dict_part)
	return NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict

def get_identity_dict(maturation_rate_dict):
	identity_dict = {}
	for (key, value) in maturation_rate_dict.items():
		for read_id in value:
			identity_dict[read_id] = key
	return identity_dict

def get_gene_usage_data(identity_dict, germline_type, NoMut_alignment_dict, Mut_alignment_dict, pic_type, germline_ids):
	outputfile1 = open("%s/%s_%s_%s_gene_usage%s.txt"%(prj_tree.data, prj_name, chain_type, germline_type, pic_type), "w")
	output_handle = csv.writer(outputfile1, delimiter="\t")
	geneusage_data, columns, max_value, gene_usage_ids = [], [], 0, []
	for germ_id in sorted(germline_ids):
		try:
			NoMut_reads_num = len(NoMut_alignment_dict[germ_id])
			NoMut_reads = NoMut_alignment_dict[germ_id]
		except KeyError:
			NoMut_reads_num = 0
			NoMut_reads = []
		try:
			Mut_reads_num = len(Mut_alignment_dict[germ_id])
			Mut_reads = Mut_alignment_dict[germ_id]
		except KeyError:
			Mut_reads_num = 0
			Mut_reads = []
		reads_num = NoMut_reads_num + Mut_reads_num
		
		Mut_reads_num_i = get_identity_list(Mut_reads, identity_dict)
		Mut_reads_num_i = continuous_subtraction_list(Mut_reads_num_i)
		Mut_reads_num_i = Mut_reads_num_i[1:] #delete < ide 75% reads
		result_line = [germ_id] + Mut_reads_num_i + [NoMut_reads_num]
		output_handle.writerow(result_line)
		geneusage_data.append(result_line[1:])
		if max_value < sum(result_line[1:]):
			max_value = sum(result_line[1:])
		columns.append(germ_id)
		gene_usage_ids.append([NoMut_reads, Mut_reads])
	outputfile1.close()
	pickle_file = '%s/%s_gene_usage_info_dump_%s_%s%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, pic_type)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = (columns, geneusage_data, gene_usage_ids)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
	outputfile1.close()
	return geneusage_data, columns
def plot_gene_usage(geneusage_data, columns, pic_type, germline_gene_list, germline_type):
	detected_genes = [x.split('*')[0] for x in columns]
	SBG = StackedBarGrapher()
	d = np.array(geneusage_data)
	total_number = np.sum(d)
	if total_number == 0:
		print "There is no %s"%germline_type
	else:
		line_numbers = np.sum(d, axis=1)
		line_number_percent_list = []
		for line_number in line_numbers:
			#line_number_percent = format((float(line_number)/float(total_number)), '.2%')
			line_number_percent = float(line_number)/float(total_number) * 100
			line_number_percent_list.append(line_number_percent)
		#print len(line_number_percent_list), len(columns)
		zip_gene_number = zip(columns, line_number_percent_list)
		#print zip_gene_number
		gene_number_dict = {}
		for item in zip_gene_number:
			gene_number_dict.setdefault(item[0].split('*')[0], []).append(item[1])
		max_len_value = 0
		for key, value in  gene_number_dict.items():
			if max_len_value < len(value):
				max_len_value = len(value)
		#print gene_number_dict
		gene_names = get_all_gene_name(germline_gene_list)
		data = np.zeros((len(gene_names),max_len_value))
		for index, gene in enumerate(gene_names):
			try:
				gene_number = gene_number_dict[gene]
				#print gene_number
				for n_index, number in enumerate(sorted(gene_number, reverse=True)):

					data[index][n_index] = number
			except KeyError:
				pass		

		d_labels = gene_names
		d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777', 
		'#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777','#2166ac', 
		'#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
		
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
			if str(t).split("'")[1] in detected_genes:
				t.set_color('blue')

		plt.savefig('%s/%s_%s_%s_gene_usage_identity%s.png'%(prj_tree.figure, prj_name, chain_type, germline_type, pic_type), dpi=300)
		del fig
		plt.close()
		
		#Plot black percentage gene usage
		black_colors = ['#000000']*len(d_colors)
		fig = plt.figure()
		ax = fig.add_subplot(111)
		SBG.stackedBarPlot(ax,
		                   data,
		                   #d_colors,
		                   black_colors,
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
			#if str(t).split("'")[1] in detected_genes:
			#	t.set_color('blue')

		plt.savefig('%s/%s_%s_%s_gene_usage_percentage%s.png'%(prj_tree.figure, prj_name, chain_type, germline_type, pic_type), dpi=300)
		del fig
		plt.close()
		#Write down gene usage percentage info
		record_gene_usage_percent(gene_number_dict, germline_type, pic_type)
		
def record_gene_usage_percent(gene_number_dict, germline_type, pic_type):
	outputfile1 = open("%s/%s_%s_%s_gene_usage_percentage%s.txt"%(prj_tree.data, prj_name, chain_type, germline_type, pic_type), "w")
	output_handle = csv.writer(outputfile1, delimiter="\t")
	for index,(key,value) in enumerate(gene_number_dict.items()):
		output_handle.writerow([key, sum(value)])
def plot_maturation_rate(maturation_rate_dict, pic_type, germline_type):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	maturation_outputfile1 = open("%s/%s_%s_%s_maturation_rate.txt"%(prj_tree.data, prj_name, chain_type, germline_type), "w")
	maturation_output_handle = csv.writer(maturation_outputfile1, delimiter="\t")
	divergence, number_reads = [], []
	for (key, value) in sorted(maturation_rate_dict.items()):
		maturation_output_handle.writerow([key, round(100-float(key), 2), len(value)])
		divergence.append(round(100-float(key), 2))
		number_reads.append(len(value))
	maturation_rate_fig = plt.figure(2)
	total_number =sum(number_reads)
	if total_number == 0:
		print "%s Maturation rate: Total number is %s"%(germline_type, total_number)
	else:
		number_reads_percent = [float(x)/total_number *100  for x in number_reads]
		plt.plot(divergence, number_reads_percent, 'bo', divergence, number_reads_percent, 'k')
		max_y_tick = math.ceil(max(number_reads_percent))
		print "Doing ytick"

		y_ticks,y_ticks_label, tick_log = [0],[0], 0
		while max_y_tick % 2 != 0:
			max_y_tick += 1
		#print "Doing ytick2"
		while tick_log <= max_y_tick:
			tick_log += 2
			y_ticks.append(tick_log*1)
			y_ticks_label.append(tick_log*1)
		#print y_ticks
		while len(y_ticks) > 14:
			if len(y_ticks) %2 == 0:
				y_ticks, y_ticks_label = y_ticks[1:][::2], y_ticks_label[1:][::2]
				y_ticks.insert(0,0), y_ticks_label.insert(0,0)
			else:
				y_ticks, y_ticks_label = y_ticks[::2], y_ticks_label[::2]
		plt.yticks(y_ticks,y_ticks_label)
		plt.xlabel('Maturation rate (%)')
		plt.ylabel('Number of reads (%)')
		plt.title("%s %s %s maturation rate"%(prj_name, chain_type, germline_type))
		plt.savefig('%s/%s_%s_%s_maturation_rate%s.png'%(prj_tree.figure, prj_name, chain_type, germline_type, pic_type))
		del fig
		plt.close()
		maturation_outputfile1.close()
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
	return sorted(set([x.split('*')[0] for x in germline_gene_list]))
def main():
	print "Begin!"
	
	#Step0 : get all reads Variable region, V gene region and CDR3 region
	origin_record_dict = SeqIO.index("%s/%s.assembled_trimed.fasta"%(prj_tree.origin, prj_name),  "fasta")
	IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
	IGBLAST_CDR3_file = "%s/%s_get_CDR3_info.txt"%(prj_tree.igblast_data, prj_name)	
	#get V-region, Variable-region, CDR3-region
	print "Trimming reads ..."
	trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file, IGBLAST_CDR3_file, origin_record_dict, chain_type)
	#Step0 over
	#Step 1:get No MUT reads and MUT reads, caculate the Ratio.
	germline_fasta = load_fasta_dict("%s/20150429-human-gl-vdj.fasta"%prj_tree.igblast_database)
	all_germline_ids = germline_fasta.keys()
	coverage80_reads_list_V, coverage80_reads_list_J = [], []
	print time.time()
	for germline_type in ('V','J'):	
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		#'''
		if germline_type == "V":
			coverage_cutoff = 5
		if germline_type == "J":
			coverage_cutoff = 0
		#os.system("rm %s/%s_get_assignment_info_dump*"%(prj_tree.tmp, prj_name))
		IGBLAST_recombanation_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
		real_reads_list, assign_result_reads_num, total_reads_num = process_recomb(IGBLAST_recombanation_file, germline_type, germline_gene_list)
		print len(real_reads_list), assign_result_reads_num, total_reads_num
		task_pool = Pool(processes = pool_size)
		IGBLAST_assignment_files = glob.glob("%s/IgBLAST_result_*_get_assignment_info.txt"%(prj_tree.igblast_data))
		for IGBLAST_assignment_file in IGBLAST_assignment_files:
			pjobs_ids = task_pool.apply_async(get_alignment_info, args=(IGBLAST_assignment_file, germline_type, real_reads_list, coverage_cutoff, ))
			#NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list, coverage_cutoff)
		print "Waiting for all subprocesses done..."
		task_pool.close()
		task_pool.join()
		#check_jobs_done(prj_name, prj_tree, "get_assignment_info", pjobs_ids)
		print 'All subprocesses done.'
		#'''
		
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict= {}, {}, {}
		get_assignment_info_dump_files = glob.glob('%s/%s_get_assignment_info_dump_%s_%s_*'%(prj_tree.tmp, prj_name, germline_type, chain_type))
		for infile in get_assignment_info_dump_files:
			NoMut_alignment_dict_part, Mut_alignment_dict_part, maturation_rate_dict_part = process_dump(infile)
			NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = combine_dict_and_list(NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict, NoMut_alignment_dict_part, Mut_alignment_dict_part, maturation_rate_dict_part)
		
		identity_dict = {}

		for (key, value) in maturation_rate_dict.items():
			for read_id in value:
				identity_dict[read_id] = key
				if germline_type == 'V':
					coverage80_reads_list_V.append(read_id)
				elif germline_type == 'J':
					coverage80_reads_list_J.append(read_id)
				else:
					pass
		print "coverage80_reads_list_V", len(coverage80_reads_list_V),len(coverage80_reads_list_J)
	#Step2 : get real reads Variable region, V gene region and CDR3 region
	print time.time()
	Variable_region_record_dict = SeqIO.index("%s/%s_%s_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
	V_gene_region_record_dict = SeqIO.index("%s/%s_%s_V_gene_region.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
	CDR3_region_record_dict = SeqIO.index("%s/%s_%s_CDR3.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
	real_reads = set(coverage80_reads_list_V) & set(coverage80_reads_list_J)
	print "real_reads", len(real_reads),  len(set(coverage80_reads_list_V)), len(set(coverage80_reads_list_J))
	real_reads_out_file_name, real_reads_out_fasta = open("%s/%s_%s_real_reads.txt"%(prj_tree.reads, prj_name, chain_type), 'w'), open("%s/%s_%s_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type), 'w')
	real_reads_out_file = csv.writer(real_reads_out_file_name, delimiter = "\t")
	V_gene_real_reads_out_fasta, CDR3_real_reads_out_fasta = open("%s/%s_%s_real_reads_V_region.fasta"%(prj_tree.reads, prj_name, chain_type), 'w'), open("%s/%s_%s_real_reads_CDR3.fasta"%(prj_tree.reads, prj_name, chain_type), 'w')
	NoCDR3_real_reads_out_fasta = open("%s/%s_%s_real_reads_No_CDR3.fasta"%(prj_tree.reads, prj_name, chain_type), 'w')
	print time.time()
	for read_id in real_reads:
		real_reads_out_file.writerow([read_id])
		SeqIO.write(Variable_region_record_dict[read_id], real_reads_out_fasta, "fasta")
		SeqIO.write(V_gene_region_record_dict[read_id], V_gene_real_reads_out_fasta, "fasta")
		try:
			SeqIO.write(CDR3_region_record_dict[read_id], CDR3_real_reads_out_fasta, "fasta")
		except KeyError:
			SeqIO.write(Variable_region_record_dict[read_id], NoCDR3_real_reads_out_fasta, "fasta")
	V_gene_real_reads_out_fasta.close()
	CDR3_real_reads_out_fasta.close()
	real_reads_out_fasta.close()
	real_reads_out_file_name.close()
	NoCDR3_real_reads_out_fasta.close()
	#Step2 Over	
	print time.time()
	for germline_type in ('V','J'):	
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		pic_type = ""
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = get_dump_data(germline_type)

		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = clear_VJcoverage80_reads(NoMut_alignment_dict, real_reads), clear_VJcoverage80_reads(Mut_alignment_dict, real_reads),  clear_VJcoverage80_reads(maturation_rate_dict, real_reads)		
		identity_dict = get_identity_dict(maturation_rate_dict)
		germline_ids = get_germ_ids(NoMut_alignment_dict.keys(), Mut_alignment_dict.keys())
		geneusage_data, columns = get_gene_usage_data(identity_dict, germline_type, NoMut_alignment_dict, Mut_alignment_dict, pic_type, germline_ids)
		#Step2: Plot gene usage and maturation rate
		plot_gene_usage(geneusage_data, columns, pic_type, germline_gene_list, germline_type)
		plot_maturation_rate(maturation_rate_dict, pic_type, germline_type)
	print time.time()
	#Step 3 Unique fasta file
	Variable_region_fasta_file = "%s/%s_%s_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type)
	Variable_region_unique_id_list, Variable_region_dict_unique = unique_fasta(Variable_region_fasta_file)
	V_region_fasta_file = "%s/%s_%s_real_reads_V_region.fasta"%(prj_tree.reads, prj_name, chain_type)
	V_region_unique_id_list, V_region_dict_unique = unique_fasta(V_region_fasta_file)
	CDR3_region_fasta_file = "%s/%s_%s_real_reads_CDR3.fasta"%(prj_tree.reads, prj_name, chain_type)
	CDR3_region_unique_id_list, CDR3_region_dict_unique = unique_fasta(CDR3_region_fasta_file)
	
	print time.time()
	Variable_region_fasta_file = SeqIO.index("%s/%s_%s_real_reads_Variable_region_unique.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
	Variable_region_unique_id_list = list(Variable_region_fasta_file.keys())
	#print Variable_region_unique_id_list
	for germline_type in ('V','J'):
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		pic_type = "_unique"
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = get_dump_data(germline_type)
		#print len(NoMut_alignment_dict)
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = clear_VJcoverage80_reads(NoMut_alignment_dict, Variable_region_unique_id_list), clear_VJcoverage80_reads(Mut_alignment_dict, Variable_region_unique_id_list), clear_VJcoverage80_reads(maturation_rate_dict, Variable_region_unique_id_list)		
		identity_dict = get_identity_dict(maturation_rate_dict)
		germline_ids = get_germ_ids(NoMut_alignment_dict.keys(), Mut_alignment_dict.keys())
		geneusage_data, columns = get_gene_usage_data(identity_dict, germline_type, NoMut_alignment_dict, Mut_alignment_dict, pic_type, germline_ids)
		#Step2: Plot gene usage and maturation rate
		plot_gene_usage(geneusage_data, columns, pic_type, germline_gene_list, germline_type)
		plot_maturation_rate(maturation_rate_dict, pic_type, germline_type)
	print time.time()
	coverage_cutoff
		
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
		
		for chain_type in ("H", "K", "L"):
			print "Processing %s chain..."%chain_type
			main()
	'''
	for chain_type in ("H", "K", "L"):
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

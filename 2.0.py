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
	reader = csv.reader(open(infile,"rU"),delimiter = "\t")
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
	return real_reads_list, count, index+1

def process_alignment(infile, germline_type, real_reads_list):
	#NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = {}, {}, [], []
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = {}, {}, {}
	print "Processing %s " %infile
	handle = open(infile,"rU")
	reader = csv.reader(handle,delimiter = "\t")

	for line in reader:
		assign_result = MyAlignment(line) #Use my_tools's MyAlignment CLASS
		coverage_rate = assign_result.coverage_rate
		if assign_result.assign_type == germline_type:
			#if total % 1000 == 0:
			#	print "%s Done!"%total
			if assign_result.identity == float(100.00) and assign_result.sstart == 1: #and assign_result.send == assign_result.slen:
				#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
				#NoMut_reads_list_all.append(assign_result.qid)
				#NoMut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#NoMut_reads_list.append(assign_result.qid)
					NoMut_alignment_dict.setdefault(assign_result.sid, []).append(assign_result.qid)
					maturation_rate_dict.setdefault(assign_result.identity, []).append(assign_result.qid)
			elif coverage_rate >= 80:
				#Mut_reads_list_all.append(assign_result.qid)
				#Mut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#Mut_reads_list.append(assign_result.qid)
					Mut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
					maturation_rate_dict.setdefault(assign_result.identity, []).append(assign_result.qid)
				
		#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)	
	#print total,len(NoMut_reads_list), len(Mut_reads_list), count
	writer.writerow(result)
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

def get_alignment_info(IGBLAST_assignment_file, germline_type, real_reads_list):
	index = IGBLAST_assignment_file.split('/')[-1].split('_')[2]
	NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list)
	pickle_file = '%s/%s_get_assignment_info_dump_%s_%s'%(prj_tree.tmp, prj_name, germline_type, index)
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
	Mut_reads_num_i = [[x, 0] for x in range(80, 101)]
	for read_id in Mut_reads:
		identity = identity_dict[read_id]
		for i in range(80,101):	
			if identity <= i:
				Mut_reads_num_i[i - 80][1] += 1
	return Mut_reads_num_i

def get_max_y_tick(max_value):
	return int(math.ceil(float(max_value)/10000))

def main():
	print "Begin!"
	origin_record_dict = SeqIO.index("%s/%s.assembled_trimed.fasta"%(prj_tree.origin, prj_name),  "fasta")

	IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
	IGBLAST_CDR3_file = "%s/%s_get_CDR3_info.txt"%(prj_tree.igblast_data, prj_name)
	
	#get V-region, Variable-region, CDR3-region
	trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file, IGBLAST_CDR3_file, origin_record_dict)

	
	
	#Step 1:get No MUT reads and MUT reads, caculate the Ratio.
	germline_fasta = load_fasta_dict("%s/IgBLAST_database/20150429-human-gl-vdj.fasta"%prj_tree.igblast_database)
	all_germline_ids = germline_fasta.keys()
	coverage80_reads_list_V, coverage80_reads_list_J = [], []
	for germline_type in ('V','J'):	
		if germline_type == 'V' and chain_type == "H":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGHV']
		if germline_type == 'J' and chain_type == "H":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGHJ']
		
		if germline_type == 'V' and chain_type == "K":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGKV']
		if germline_type == 'J' and chain_type == "K":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGKJ']
		
		if germline_type == 'V' and chain_type == "L":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGLV']
		if germline_type == 'J' and chain_type == "L":
			germline_gene_list = HUMAN_GERMLINE['HUMANIGLJ']
		'''
		#os.system("rm %s/%s_get_assignment_info_dump*"%(prj_tree.tmp, prj_name))
		IGBLAST_recombanation_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
		real_reads_list, assign_result_reads_num, total_reads_num = process_recomb(IGBLAST_recombanation_file, germline_type, germline_gene_list)
		print len(real_reads_list), assign_result_reads_num, total_reads_num
		#sys.exit(0)
		task_pool = Pool(processes = pool_size)
		IGBLAST_assignment_files = glob.glob("%s/IgBLAST_result_*_get_assignment_info.txt"%(prj_tree.igblast_data))
		for IGBLAST_assignment_file in IGBLAST_assignment_files:
			pjobs_ids = task_pool.apply_async(get_alignment_info, args=(IGBLAST_assignment_file, germline_type, real_reads_list,))
			#NoMut_alignment_dict, Mut_alignment_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list)
		print "Waiting for all subprocesses done..."
		task_pool.close()
		task_pool.join()
		#check_jobs_done(prj_name, prj_tree, "get_assignment_info", pjobs_ids)
		print 'All subprocesses done.'
		'''
		
		NoMut_alignment_dict, Mut_alignment_dict, maturation_rate_dict= {}, {}, {}
		get_assignment_info_dump_files = glob.glob('%s/%s_get_assignment_info_dump_%s_*'%(prj_tree.tmp, prj_name, germline_type))
		for infile in get_assignment_info_dump_files:
			print infile
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
		'''
		germline_ids = get_germ_ids(NoMut_alignment_dict.keys(), Mut_alignment_dict.keys())
		outputfile1 = "%s/%s_%s_gene_usage.txt"%(prj_tree.data, prj_name, germline_type)
		output_handle = csv.writer(open(outputfile1, "w"), delimiter="\t")
		outputfile2 = "%s/%s_%s_gene_usage_sorted.txt"%(prj_tree.data, prj_name, germline_type)
		output_handle_sorted = csv.writer(open(outputfile2, "w"), delimiter="\t")
		geneusage_data, columns, max_value = [], [], 0
		for germ_id in sorted(germline_ids):
			try:
				NoMut_reads_num = len(NoMut_alignment_dict[germ_id])
			except KeyError:
				NoMut_reads_num = 0
			try:
				Mut_reads_num = len(Mut_alignment_dict[germ_id])
				Mut_reads = Mut_alignment_dict[germ_id]
			except KeyError:
				Mut_reads_num = 0
				Mut_reads = []
			reads_num = NoMut_reads_num + Mut_reads_num
			Mut_reads_num_i = get_identity_list(Mut_reads, identity_dict)
			Mut_reads_num_i = continuous_subtraction_list(Mut_reads_num_i)
			result_line = [germ_id] + Mut_reads_num_i + [NoMut_reads_num]
			output_handle.writerow(result_line)
			geneusage_data.append(result_line[1:])
			if max_value < sum(result_line[1:]):
				max_value = sum(result_line[1:])
			columns.append(germ_id)
			if germ_id.split('*')[0] in germline_gene_list:
				output_handle_sorted.writerow([germ_id] + Mut_reads_num_i + [NoMut_reads_num])
		
		#Step2: Plot gene usage and maturation rate

		SBG = StackedBarGrapher()
		d = np.array(geneusage_data)
		d_labels = columns
		d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
		rows = ['%d%%' % x for x in range(80, 101)]
		rows.insert(0, "<80%")
		colors = plt.cm.BuPu(np.linspace(0, 0.5, len(rows)))
		d_colors = colors
		fig = plt.figure(1)
		max_y_tick = get_max_y_tick(max_value)
		y_ticks, tick_log = [[0],[0]], 1
		while tick_log <= max_y_tick:
			y_ticks[0].append(tick_log*10000)
			y_ticks[1].append(tick_log*10)
			tick_log += 1
		ax = fig.add_subplot(111)
		SBG.stackedBarPlot(ax,
		                   d,
		                   d_colors,
		                   xLabels=d_labels,
		                   yTicks= y_ticks,
		                   scale=False,
						   ylabel = 'Number of reads (*1000)',
						   gap =.2
		                  )
		plt.title("%s_%s_gene_usage"%(prj_name, germline_type))

		fig.subplots_adjust(bottom=0.4)
		plt.tight_layout()
		
		plt.savefig('%s/%s_%s_gene_usage.png'%(prj_tree.figure, prj_name, germline_type))
		maturation_outputfile1 = "%s/%s_%s_maturation_rate.txt"%(prj_tree.data, prj_name, germline_type)
		maturation_output_handle = csv.writer(open(maturation_outputfile1, "w"), delimiter="\t")
		maturation_outputfile2 = "%s/%s_%s_maturation_rate.txt"%(prj_tree.data, prj_name, germline_type)
		maturation_output_handle_sorted = csv.writer(open(maturation_outputfile2, "w"), delimiter="\t")
		divergence, number_reads = [], []
		for (key, value) in sorted(maturation_rate_dict.items()):
			maturation_output_handle.writerow([key, round(100-float(key), 2), len(value)])
			divergence.append(round(100-float(key), 2))
			number_reads.append(len(value))
		maturation_rate_fig = plt.figure(2)
		plt.plot(divergence, number_reads, 'bo', divergence, number_reads, 'k')
		max_y_tick = get_max_y_tick(max(number_reads))
		
		#if max_y_tick % 2 != 0:
		#	max_y_tick += 1
		y_ticks, tick_log = [[0],[0]], 0
		while tick_log <= max_y_tick:
			tick_log += 4
			y_ticks[0].append(tick_log*10000)
			y_ticks[1].append(tick_log*10)
			
		plt.yticks(y_ticks[0],y_ticks[1])
		plt.xlabel('Maturation rate (%)')
		plt.ylabel('Number of reads (*1000)')
		plt.title("%s_%s_maturation_rate"%(prj_name, germline_type))
		plt.savefig('%s/%s_%s_maturation_rate.png'%(prj_tree.figure, prj_name, germline_type))
		'''
	Variable_region_record_dict = SeqIO.index("%s/%s_Variable_region.fasta"%(prj_tree.reads, prj_name), "fasta")
	V_gene_region_record_dict = SeqIO.index("%s/%s_V_gene_region.fasta"%(prj_tree.reads, prj_name), "fasta")
	CDR3_region_record_dict = SeqIO.index("%s/%s_CDR3.fasta"%(prj_tree.reads, prj_name), "fasta")
	real_reads = set(coverage80_reads_list_V) & set(coverage80_reads_list_J)
	real_reads_out_file, real_reads_out_fasta = csv.writer(open("%s/%s_real_reads.txt"%(prj_tree.reads, prj_name), 'w'), delimiter = "\t"), open("%s/%s_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name), 'w')
	V_gene_real_reads_out_fasta, CDR3_real_reads_out_fasta = open("%s/%s_real_reads_V_region.fasta"%(prj_tree.reads, prj_name), 'w'), open("%s/%s_real_reads_CDR3.fasta"%(prj_tree.reads, prj_name), 'w')
	NoCDR3_real_reads_out_fasta = open("%s/%s_real_reads_No_CDR3.fasta"%(prj_tree.reads, prj_name), 'w')
	for read_id in real_reads:
		real_reads_out_file.writerow([read_id])
		SeqIO.write(Variable_region_record_dict[read_id], real_reads_out_fasta, "fasta")
		SeqIO.write(V_gene_region_record_dict[read_id], V_gene_real_reads_out_fasta, "fasta")
		try:
			SeqIO.write(CDR3_region_record_dict[read_id], CDR3_real_reads_out_fasta, "fasta")
		except KeyError:
			SeqIO.write(Variable_region_record_dict[read_id], NoCDR3_real_reads_out_fasta, "fasta")
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	if "K" in prj_name:
		chain_type = "K"
	elif "L" in prj_name:
		chain_type = "L"
	else:
		chain_type = "H"
	main()
	print "Finished"

#!/usr/bin/env python
# encoding: utf-8
"""
2.0-process-igblast-output.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
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

def process_recomb(infile, germline_type):
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
		if recomb_result.productive == "Yes" and recomb_result.cover_vj == "Yes":
			real_reads_list.append(recomb_result.qid)
	return real_reads_list, count, index+1
def process_alignment(infile, germline_type, real_reads_list):
	#NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = {}, {}, [], []
	NoMut_alignment_dict, Mut_alignment_dict = {}, {}
	print "Processing %s " %infile
	reader = csv.reader(open(infile,"rU"),delimiter = "\t")
	total, count = 0, 0
	for line in reader:
		assign_result = MyAlignment(line) #Use my_tools's MyAlignment CLASS
		coverage_rate = assign_result.coverage_rate
		if assign_result.assign_type == germline_type:
			total += 1
			#if total % 1000 == 0:
			#	print "%s Done!"%total
			if assign_result.identity == float(100.00) and assign_result.sstart == 1 and assign_result.send == assign_result.slen:
				#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
				#NoMut_reads_list_all.append(assign_result.qid)
				#NoMut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#NoMut_reads_list.append(assign_result.qid)
					NoMut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
					if None == NoMut_alignment_dict[assign_result.sid]:
						print NoMut_alignment_dict
						sys.exit(0)
			elif coverage_rate >= 80:
				#Mut_reads_list_all.append(assign_result.qid)
				#Mut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					#Mut_reads_list.append(assign_result.qid)
					Mut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
			else:
				count += 1
				
		#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
		
	#print total,len(NoMut_reads_list), len(Mut_reads_list), count
	return NoMut_alignment_dict, Mut_alignment_dict

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
	NoMut_alignment_dict, Mut_alignment_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list)
	pickle_file = '%s/%s_get_assignment_info_dump_%s'%(prj_tree.tmp, prj_name, index)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = (NoMut_alignment_dict, Mut_alignment_dict)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()

def process_dump(infile):
	f = open(infile, 'rb')
	pickle_tuple = pickle.load(f)
	#NoMut_alignment_dict, Mut_alignment_dict = pickle_tuple[0], pickle_tuple[1]
	NoMut_alignment_dict, Mut_alignment_dict, NoMut_reads_list, Mut_reads_list, NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = pickle_tuple[0], pickle_tuple[1], pickle_tuple[2], pickle_tuple[3], pickle_tuple[4], pickle_tuple[5], pickle_tuple[6], pickle_tuple[7]
	#for (key, value) in NoMut_alignment_dict.items():
		#if value == None:
		#print NoMut_alignment_dict
		#sys.exit(0)
	f.close()
	return	NoMut_alignment_dict, Mut_alignment_dict
def combine_dict_and_list(NoMut_alignment_dict, Mut_alignment_dict, NoMut_alignment_dict_part, Mut_alignment_dict_part):
	NoMut_alignment_dict_combined, Mut_alignment_dict_combined = {}, {}
	germline_ids = get_germ_ids(NoMut_alignment_dict.keys(), NoMut_alignment_dict_part.keys())
	for germline_id in germline_ids:
		try:
			NoMut_reads = NoMut_alignment_dict[germline_id]
		except KeyError:
			NoMut_reads = []
		try:
			NoMut_reads_part = NoMut_alignment_dict_part[germline_id]
		except KeyError:
			NoMut_reads_part = []
		NoMut_alignment_dict_combined[germline_id] = NoMut_reads + NoMut_reads_part
	
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
	return NoMut_alignment_dict_combined, Mut_alignment_dict_combined
def main():
	print "Begin!"
	task_pool = Pool(processes = pool_size-1)
	'''
	IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
	trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file)
	'''
	#Step 1:get No MUT reads and MUT reads, caculate the Ratio.
	germline_fasta = load_fasta_dict("%s/IgBLAST_database/20150429-human-gl-vdj.fasta"%prj_tree.igblast_database)
	all_germline_ids = germline_fasta.keys()
	for germline_type in ('V','J'):
		outputfile = "%s/%s_%s_gene_usage.txt"%(prj_tree.data, prj_name, germline_type)
		output_handle = csv.writer(open(outputfile, "w"), delimiter="\t")
		'''
		IGBLAST_recombanation_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
		real_reads_list, assign_result_reads_num, total_reads_num = process_recomb(IGBLAST_recombanation_file, germline_type)
		print len(real_reads_list), assign_result_reads_num, total_reads_num
		#sys.exit(0)
		
		IGBLAST_assignment_files = glob.glob("%s/IgBLAST_result_*_get_assignment_info.txt"%(prj_tree.igblast_data))
		for IGBLAST_assignment_file in IGBLAST_assignment_files:
			#pjobs_ids = task_pool.apply_async(get_alignment_info, args=(IGBLAST_assignment_file, germline_type, real_reads_list,))
			NoMut_alignment_dict, Mut_alignment_dict = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list)
		print "Waiting for all subprocesses done..."
		task_pool.close()
		task_pool.join()
		#check_jobs_done(prj_name, prj_tree, "get_assignment_info", pjobs_ids)
		print 'All subprocesses done.'
		'''
		NoMut_alignment_dict, Mut_alignment_dict= {}, {}
		get_assignment_info_dump_files = glob.glob('%s/%s_get_assignment_info_dump_*'%(prj_tree.tmp, prj_name))
		for infile in get_assignment_info_dump_files:
			print infile
			NoMut_alignment_dict_part, Mut_alignment_dict_part = process_dump(infile)
			NoMut_alignment_dict, Mut_alignment_dict = combine_dict_and_list(NoMut_alignment_dict, Mut_alignment_dict, NoMut_alignment_dict_part, Mut_alignment_dict_part)
		print "Mark3", [(key, len(value)) for (key, value) in sorted(NoMut_alignment_dict.items())]
		print "Mark4", [(key, len(value)) for (key, value) in sorted(Mut_alignment_dict.items())]
		'''
		print len(NoMut_alignment_dict.keys()), len(NoMut_reads_list)
		print len(Mut_alignment_dict.keys()), len(Mut_reads_list)
		print len(NoMut_alignment_dict_all.keys()), len(NoMut_reads_list_all)
		print len(Mut_alignment_dict_all.keys()), len(Mut_reads_list_all)
		sys.exit(0)
		'''
		germline_ids = get_germ_ids(NoMut_alignment_dict.keys(), Mut_alignment_dict.keys())
		for germ_id in germline_ids:
			try:
				NoMut_reads_num = len(NoMut_alignment_dict[germ_id])
			except KeyError:
				NoMut_reads_num = 0
			try:
				Mut_reads_num = len(Mut_alignment_dict[germ_id])
			except KeyError:
				Mut_reads_num = 0
			reads_num = NoMut_reads_num + Mut_reads_num
			output_handle.writerow([germ_id, NoMut_reads_num, Mut_reads_num, reads_num])
		print len(real_reads_list), assign_result_reads_num, total_reads_num
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	
	main()
	print "Finished"

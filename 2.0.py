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
	NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = {}, {}, [], []
	NoMut_alignment_dict, Mut_alignment_dict, NoMut_reads_list, Mut_reads_list = {}, {}, [], []
	print "Processing %s " %infile
	reader = csv.reader(open(infile,"rU"),delimiter = "\t")
	total, count = 0, 0
	for line in reader:
		assign_result = MyAlignment(line) #Use my_tools's MyAlignment CLASS
		coverage_rate = assign_result.coverage_rate
		if assign_result.assign_type == germline_type:
			total += 1
			if total % 1000 == 0:
				print "%s Done!"%total
			if assign_result.identity == float(100.00)  and assign_result.assign_type == germline_type and assign_result.sstart == 1 and assign_result.send == assign_result.slen:
				#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
				NoMut_reads_list_all.append(assign_result.qid)
				NoMut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					NoMut_reads_list.append(assign_result.qid)
					NoMut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
			elif coverage_rate >= 80 and assign_result.assign_type == germline_type:
				Mut_reads_list_all.append(assign_result.qid)
				Mut_alignment_dict_all.setdefault(assign_result.sid,[]).append(assign_result.qid)
				if assign_result.qid in real_reads_list:
					Mut_reads_list.append(assign_result.qid)
					Mut_alignment_dict.setdefault(assign_result.sid,[]).append(assign_result.qid)
			else:
				count += 1
				
		#print assign_result.qid,assign_result.sid,assign_result.strand, assign_result.qseq, list(assign_result.qseq)
		
	#print total,len(NoMut_reads_list), len(Mut_reads_list), count
	return NoMut_alignment_dict, Mut_alignment_dict, NoMut_reads_list, Mut_reads_list, NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all	


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

def main():
	print "Begin!"
	
	#Step 1:get No MUT reads and MUT reads, caculate the Ratio.
	
	for germline_type in ('V','J'):
		'''
		IGBLAST_recombanation_file = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
		real_reads_list, assign_result_reads_num, total_reads_num = process_recomb(IGBLAST_recombanation_file, germline_type)
		print len(real_reads_list), assign_result_reads_num, total_reads_num
		#sys.exit(0)
		'''
		IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
		trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file)
		#IGBLAST_assignment_files = glob.glob("%s/IgBLAST_result_*_get_assignment_info.txt"%(prj_tree.igblast_data))
		#NoMut_alignment_dict, Mut_alignment_dict, NoMut_reads_list, Mut_reads_list, NoMut_alignment_dict_all, Mut_alignment_dict_all, NoMut_reads_list_all, Mut_reads_list_all = process_alignment(IGBLAST_assignment_file, germline_type, real_reads_list)
		'''
		print len(NoMut_alignment_dict.keys()), len(NoMut_reads_list)
		print len(Mut_alignment_dict.keys()), len(Mut_reads_list)
		print len(NoMut_alignment_dict_all.keys()), len(NoMut_reads_list_all)
		print len(Mut_alignment_dict_all.keys()), len(Mut_reads_list_all)
		sys.exit(0)
		'''
		value_list = def_score_function(alignment_dict)
		
		NoMut_line_writer = csv.writer(open("NoMut_line_germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		NoMut_line_writer.writerows(NoMut_result_list)
		
		unique_ID = calculate_score(value_list)
		writer = csv.writer(open("germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		
		NoMut_value_list = def_score_function(NoMut_alignment_dict)
		NoMut_unique_ID = calculate_score(NoMut_value_list)
		NoMut_writer = csv.writer(open("NoMut_germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		
		for item in sorted(unique_ID):
			writer.writerow(item)
		
		for item in sorted(NoMut_unique_ID):
			NoMut_writer.writerow(item)

if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	task_pool = Pool(processes = pool_size-1)
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	
	main()
	print "Finished"

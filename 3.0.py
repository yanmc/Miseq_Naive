#!/usr/bin/env python
# encoding: utf-8
"""
3.0.py 

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
		if len(value) > 1:
			sorted_value = sorted(value, key = lambda z : z[1], reverse=True)
			if float(sorted_value[1][1])/float(sorted_value[0][1]) < 0.001:
				sorted_value = [sorted_value[0][0]]
			else:
				sorted_value = [x[0] for x in sorted_value[:2]]
		else:
			sorted_value = [value[0][0]]
		result_dict[key] = sorted_value
	result_file_name = open('%s/%s_%s_%s_right_allele_usage%s.txt'%(prj_tree.data, prj_name, chain_type, germline_type, pic_type),"w")
	result_file = csv.writer(result_file_name, delimiter="\t")
	result_file.writerow(["Germline gene", "Allele 1", "Allele 2"])
	max_freq_allele_dict = {}
	for germline_gene in sorted(germline_gene_list):
		try:
			value = result_dict[germline_gene]
			if len(value) == 1:
				result_file.writerow([germline_gene, value[0].split('*')[1], 0])
				max_freq_allele_dict[germline_gene] = value[0]
			else:
				result_file.writerow([germline_gene, value[0].split('*')[1], value[1].split('*')[1]])
				max_freq_allele_dict[germline_gene] = value[0]
		except:
			pass
	result_file_name.close()
	return max_freq_allele_dict
def get_novel_allele(assignment_dict, result_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta):
	gene_reads_id_dict = {}
	for (key, value) in result_ids_dict.items():
		gene_reads_id_dict.setdefault(key.split('*')[0], []).extend(value)
	for key,value in gene_reads_id_dict.items():
		print key, type(value), len(value)

	for genmline_gene in germline_gene_list:
		print "Processing gene: %s"%genmline_gene
		try:
			max_freq_allele = max_freq_allele_dict[genmline_gene]
		except KeyError:
			continue
		all_reads_ids_list = gene_reads_id_dict[genmline_gene]
		get_mutation_patterns_v2(assignment_dict, all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta, germline_type)
		#get_mutation_patterns( all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta)
		#get_mutation_spectrum(all_reads_ids_list, max_freq_allele)

def get_mutation_patterns(all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta):
	mutation_patterns_dict = {}
	print type(all_reads_ids_list), len(all_reads_ids_list)
	for read_id in all_reads_ids_list:
		ref_seq_id, test_seq_id = max_freq_allele, read_id
		#print ref_seq_id, type(test_seq_id),len(test_seq_id), type(germline_fasta), type(unique_real_reads_fasta)
		ref_seqrecord = germline_fasta[ref_seq_id]
		test_seqrecord = unique_real_reads_fasta[test_seq_id]
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
	for (key, value) in mutation_patterns_dict.items():
		mutation_patterns_group.setdefault(value[0], []).append((key, value[1]))
	data = np.zeros((len(mutation_patterns_group)-1, len(germline_fasta[max_freq_allele])))
	for index, (group_number, value) in enumerate(mutation_patterns_group.items()):
		if index != 0:
			position_mut_list = [item[1]  for item in value]
			for position_list in position_mut_list:
				for position_record in position_list:
					position = position_record[0]
					data[index-1][position-1] += 1
	ref_seq_id_name = max_freq_allele.split('*')[0].replace('/','')
	mutation_patterns_file = open('%s/%s_%s_mutation_patterns.txt'%(prj_tree.data, prj_name, ref_seq_id_name), 'w')
	mutation_patterns_writer = csv.writer(mutation_patterns_file, delimiter = "\t")
	#for line in data:
	#	mutation_patterns_writer.writerow(line)
	mutation_patterns_writer.writerows(data)
	mutation_patterns_file.close()
	return data
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
	mutation_patterns_result = []
	for index, (x, y) in enumerate(zip_seqs):
		if x != y:
			position = index + 1
			mutation_patterns_result.append((position, x, y))
	return mutation_patterns_result, mismatches
def caculate_mutation_patterns(clustalw_result, read_id, mutation_patterns_dict):
	result_file = open('%s_mutation_patterns.txt'%prj_name,'a+')
	result = csv.writer(result_file,delimiter = '\t')
	fs = glob.glob(clustalw_result)
	for infile in fs:
		#print "processing %s"%infile
		mutation_patterns_result, mismatches = parse_pair_clustal2(infile)
		mutation_patterns_dict[read_id] = (mismatches, mutation_patterns_result)
	result_file.close()
	return mutation_patterns_dict
def do_clustalw(file_for_clustalw):
	infiles = glob.glob(file_for_clustalw)

	#clustalw_exe = r"/Applications/clustalw-2.1-macosx/clustalw2"
	clustalw_exe = r"/zzh_gpfs/apps/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
	assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
	for in_file in infiles:
		#print "Processing %s ......."%in_file
		in_file = in_file.replace('&','\&')
		in_file = in_file.replace('*','\*')
		clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file)
		stdout, stderr = clustalw_cline()

def chunks(arr, n):
	return [arr[i:i + n] for i in range(0, len(arr), n)]

def get_mutation_patterns_dict(assignment_dict, all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta, germline_type,file_index):
	mutation_patterns_dict = {}
	for read_index, read_id in enumerate(all_reads_ids_list):
		if read_index % 10000 == 0:
			print "%s has been processed..."%read_index
		ref_seq_id, test_seq_id = max_freq_allele, read_id
		assignment_result_VDJ = assignment_dict[test_seq_id]
		for item in assignment_result_VDJ:
			if item[0] == germline_type and item[4] == ref_seq_id:
				assignment_result = item
				assign_start = int(assignment_result[1])
				query_seq, ref_seq = assignment_result[2], assignment_result[3]
				ref_seq, query_seq, trimmed_ref_len = remove_ref_insertion(ref_seq, query_seq)
				ref_query_seq_zip = zip(ref_seq, query_seq)
				mismatches	= sum([x != y for x, y in ref_query_seq_zip])
				mutation_patterns_result = []
				for index, (x, y) in enumerate(ref_query_seq_zip):
					if x != y:
						position = index
						mutation_patterns_result.append((assign_start+position, x, y))
				mutation_patterns_dict[read_id] = (mismatches, mutation_patterns_result)
		if read_id in mutation_patterns_dict.keys():
			pass
		else:
			#print ref_seq_id, type(test_seq_id),len(test_seq_id), type(germline_fasta), type(unique_real_reads_fasta)
			ref_seqrecord = germline_fasta[ref_seq_id]
			test_seqrecord = unique_real_reads_fasta[test_seq_id]
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
	pickle_file = '%s/%s_get_mutation_patterns_dict_dump_%s_%s_%s_%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, max_freq_allele.split('*')[0],file_index)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = mutation_patterns_dict
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
def get_mutation_patterns_v2(assignment_dict, all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta, germline_type):
	
	print type(all_reads_ids_list), len(all_reads_ids_list)
	all_reads_ids_list_chunks = chunks(all_reads_ids_list, 500)
	task_pool = Pool(processes = pool_size)
	
	for (index, all_reads_ids_list) in enumerate(all_reads_ids_list_chunks):
		print "Processing %s part"%index
		#get_mutation_patterns_dict(assignment_dict, all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta, germline_type,index)
		pjobs_ids = task_pool.apply_async(get_mutation_patterns_dict, args=(assignment_dict, all_reads_ids_list, max_freq_allele,  germline_fasta, unique_real_reads_fasta, germline_type,index,))
	print "Waiting for all subprocesses done..."
	task_pool.close()
	task_pool.join()
	#check_jobs_done(prj_name, prj_tree, "get_assignment_info", pjobs_ids)
	print 'All subprocesses done.'
	#'''
	mutation_patterns_dict = {}		
	pickle_files = glob.glob('%s/%s_get_mutation_patterns_dict_dump_%s_%s_%s_*'%(prj_tree.tmp, prj_name, germline_type, chain_type, max_freq_allele.split('*')[0]))
	for pickle_file in pickle_files:
		mutation_patterns_dict_part = process_dump(pickle_file)
		mutation_patterns_dict = combine_dict_and_list(mutation_patterns_dict, mutation_patterns_dict_part)
		
	mutation_patterns_group = {}
	for (key, value) in mutation_patterns_dict.items():
		mutation_patterns_group.setdefault(value[0], []).append((key, value[1]))
		
	
	ref_seq_id_name = max_freq_allele.split('*')[0].replace('/','')
	pickle_file = '%s/%s_mutation_pattrens_dump_%s_%s_%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, ref_seq_id_name)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = (mutation_patterns_group)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
	
	data = np.zeros((len(mutation_patterns_group)-1, len(germline_fasta[max_freq_allele])))
	for index, (group_number, value) in enumerate(mutation_patterns_group.items()):
		if index != 0:
			position_mut_list = [item[1]  for item in value]
			for position_list in position_mut_list:
				for position_record in position_list:
					position = position_record[0]
					try:
						data[index-1][position-1] += 1
					except IndexError:
						print index, position, group_number, value
	mutation_patterns_file = open('%s/%s_%s_mutation_patterns.txt'%(prj_tree.data, prj_name, ref_seq_id_name), 'w')
	mutation_patterns_writer = csv.writer(mutation_patterns_file, delimiter = "\t")
	#for line in data:
	#	mutation_patterns_writer.writerow(line)
	mutation_patterns_writer.writerows(data)
	mutation_patterns_file.close()
	return data	

def process_dump(infile):
	f = open(infile, 'rb')
	pickle_tuple = pickle.load(f)
	mutation_patterns_dict = pickle_tuple
	f.close()
	return	mutation_patterns_dict

def get_germ_ids(x, y):
	return set(x)|set(y)
def combine_dict_and_list(NoMut_alignment_dict, NoMut_alignment_dict_part):
	
	NoMut_alignment_dict_combined = combine_dict_part(NoMut_alignment_dict, NoMut_alignment_dict_part)
	return NoMut_alignment_dict_combined

def combine_dict_part(dict1, dict2):
	dictMerged2=dict(dict1, **dict2)
	
	return dictMerged2
def load_assignment_dict(IGBLAST_assignment_file):
	assignment_dict = {}
	infile = open(IGBLAST_assignment_file,"rU")
	reader = csv.reader(infile,delimiter = "\t")
	for index, line in enumerate(reader):
		assign_result = MyAlignment(line)
		assignment_dict.setdefault(assign_result.qid, []).append([assign_result.assign_type, assign_result.sstart, assign_result.qseq, assign_result.sseq, assign_result.sid])
	infile.close()
	return assignment_dict
def main():
	#germline_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta")
	#unique_real_reads_fasta = load_fasta_dict("%s/%s_H_real_reads_V_region_unique.fasta"%(prj_tree.reads, prj_name))
	germline_fasta = SeqIO.to_dict(SeqIO.parse("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta", "fasta"))
	unique_real_reads_fasta = SeqIO.to_dict(SeqIO.parse("%s/%s_%s_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta"))
	
	IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
	assignment_dict = load_assignment_dict(IGBLAST_assignment_file)
	
	manager = Manager()
	assignment_dict = manager.dict(assignment_dict)
	germline_fasta = manager.dict(germline_fasta)
	unique_real_reads_fasta = manager.dict(unique_real_reads_fasta)
	
	for germline_type in ('V','J'):
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		pickle_file = '%s/%s_gene_usage_info_dump_%s_%s%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, pic_type)
		f_IgH = open(pickle_file, 'rb')
		pickle_tuple_IgH = pickle.load(f_IgH)
		columns_IgH, geneusage_data_IgH, gene_usage_ids_IgH = pickle_tuple_IgH[0], pickle_tuple_IgH[1], pickle_tuple_IgH[2]
		f_IgH.close()
		print type(gene_usage_ids_IgH), len(gene_usage_ids_IgH)
		IgH_geneusage_dict, IgH_geneusage_ids_dict = {}, {}
		for i in range(len( pickle_tuple_IgH[0])):
			print "detect length:", type(gene_usage_ids_IgH[i][0]), type(gene_usage_ids_IgH[i][1]), len(gene_usage_ids_IgH[i][0]), len(gene_usage_ids_IgH[i][1])
			IgH_geneusage_dict[columns_IgH[i]] = sum(geneusage_data_IgH[i])
			IgH_geneusage_ids_dict[columns_IgH[i]] = set(gene_usage_ids_IgH[i][0]) | set(gene_usage_ids_IgH[i][1])
			print "detect length:", len(gene_usage_ids_IgH[i][0]), len(gene_usage_ids_IgH[i][1]), len(IgH_geneusage_ids_dict[columns_IgH[i]])
		max_freq_allele_dict = get_naive_right_allele(IgH_geneusage_dict, pic_type, germline_gene_list, germline_type)
		print "get_novel_allele"
		get_novel_allele(assignment_dict, IgH_geneusage_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta)
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()
	pic_type = ""#"_unique"
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
			main()
	end = time.time()
	print prj_name, end-start
	print "Finished"

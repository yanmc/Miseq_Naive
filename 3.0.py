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
	result_file = csv.writer(open('%s/%s_%s_%s_right_allele_usage%s.txt'%(prj_tree.data, prj_name, chain_type, germline_type, pic_type),"w"), delimiter="\t")
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
	return max_freq_allele_dict
def get_novel_allele(result_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta):
	gene_reads_id_dict = {}
	for (key, value) in result_ids_dict.items():
		gene_reads_id_dict.setdefault(key.split('*')[0], []).extend(value)
	for key,value in gene_reads_id_dict.items():
		print key, type(value), len(value)

	for genmline_gene in germline_gene_list:
		try:
			max_freq_allele = max_freq_allele_dict[genmline_gene]
		except KeyError:
			continue
		all_reads_ids_list = gene_reads_id_dict[genmline_gene]
		get_mutation_patterns(all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta)
		#get_mutation_spectrum(all_reads_ids_list, max_freq_allele)

def get_mutation_patterns(all_reads_ids_list, max_freq_allele, germline_fasta, unique_real_reads_fasta):
	mutation_patterns_dict = {}
	print type(all_reads_ids_list), len(all_reads_ids_list)
	for read_id in all_reads_ids_list:
		ref_seq_id, test_seq_id = max_freq_allele, read_id
		print ref_seq_id, type(test_seq_id),len(test_seq_id), type(germline_fasta), type(unique_real_reads_fasta)
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
		print group_number, value

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
	#germline_fasta = load_fasta_dict("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta")
	#unique_real_reads_fasta = load_fasta_dict("%s/%s_H_real_reads_V_region_unique.fasta"%(prj_tree.reads, prj_name))
	germline_fasta = SeqIO.index("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta", "fasta")
	unique_real_reads_fasta = SeqIO.index("%s/%s_%s_real_reads_V_region.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
	for germline_type in ('V','J'):
		germline_gene_list = get_germline_gene(germline_type, chain_type)
		pickle_file = '%s/%s_gene_usage_info_dump_%s_%s%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, pic_type)
		f_IgH = open(pickle_file, 'rb')
		pickle_tuple_IgH = pickle.load(f_IgH)
		columns_IgH, geneusage_data_IgH, gene_usage_ids_IgH = pickle_tuple_IgH[0], pickle_tuple_IgH[1], pickle_tuple_IgH[2]
		print type(gene_usage_ids_IgH), len(gene_usage_ids_IgH)
		IgH_geneusage_dict, IgH_geneusage_ids_dict = {}, {}
		for i in range(len( pickle_tuple_IgH[0])):
			print "detect length:", type(gene_usage_ids_IgH[i][0]), type(gene_usage_ids_IgH[i][1]), len(gene_usage_ids_IgH[i][0]), len(gene_usage_ids_IgH[i][1])
			IgH_geneusage_dict[columns_IgH[i]] = sum(geneusage_data_IgH[i])
			IgH_geneusage_ids_dict[columns_IgH[i]] = set(gene_usage_ids_IgH[i][0]) | set(gene_usage_ids_IgH[i][1])
			print "detect length:", len(gene_usage_ids_IgH[i][0]), len(gene_usage_ids_IgH[i][1]), len(IgH_geneusage_ids_dict[columns_IgH[i]])
		max_freq_allele_dict = get_naive_right_allele(IgH_geneusage_dict, pic_type, germline_gene_list, germline_type)
		print "get_novel_allele"
		get_novel_allele(IgH_geneusage_ids_dict, pic_type, germline_gene_list, germline_type, max_freq_allele_dict, germline_fasta, unique_real_reads_fasta)
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()
	pic_type = "_unique"
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

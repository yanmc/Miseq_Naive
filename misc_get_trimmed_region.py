#!/usr/bin/env python
# encoding: utf-8
"""
6.1-get-trimmed-region.py -i infile -ref reference_file -org orgnism

Created by Mingchen on 2014-12-11.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  original fasta file
orgnism:  human,mouse or rabbit
referecce_file:  the file contain variable region start and end info eg: human_get_vdj.txt 

note: already trans the reversed reads 
"""
import glob, csv, Bio, re,os,sys
from Bio import SeqIO
from mytools import *
from Bio import SeqIO
def tanslate_file(infile):
	outfile = open("%s_AA.fasta"%infile.split(".")[0], "w")
	handlefile = SeqIO.parse(infile,"fasta")
	for (index,record) in enumerate(handlefile):
		#print record
		aa = record.seq.translate()
		count, trim_seq = 0, record.seq
		while "*" in aa and count < 3:
			#print aa
			trim_seq = trim_seq[1:]
			aa = trim_seq.translate()
			count += 1
		if "*" in aa or aa == "":
			pass
		else:
			#print aa
			seqrecord = SeqRecord(aa, id = record.id, description = "Trim:%s"%count)
			SeqIO.write(seqrecord, outfile, "fasta")
def trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file, IGBLAST_CDR3_file, origin_record_dict, chain_type):
	leader_region_outfile = open("%s/%s_%s_leader_region.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	C_region_outfile = open("%s/%s_%s_C_region.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	CDR3_outfile = open("%s/%s_%s_CDR3.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	CDR3_reader = csv.reader(open(IGBLAST_CDR3_file,"rU"), delimiter = "\t")
	CDR3_start_dict = {}
	for line in CDR3_reader:
		try:
			query_id = line[0].strip().split(' ')[0]
			CDR3_start = line[2]
		except IndexError:
			continue
		CDR3_start_dict[query_id] = CDR3_start
	
	outfile = open("%s/%s_%s_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	outfile2 = open("%s/%s_%s_V_gene_region.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	outfile3 = open("%s/%s_%s_J_gene_region.fasta"%(prj_tree.reads, prj_name, chain_type),"w")
	reader = csv.reader(open(IGBLAST_assignment_file,"rU"), delimiter = "\t")
	result = ['Query_ID','query_seq','query_length','Variable_region', 'V_gene_region', 'Leader_region', 'C_region','J_region']
	assign_position_dict = {}
	for line in reader:
		assign_result = MyAlignment(line)
		if assign_result.assign_type == "V":
			if assign_result.strand	== "+":
				origin_seq = origin_record_dict[assign_result.qid]
			elif assign_result.strand	== "-":
				origin_seq = origin_record_dict[assign_result.qid].reverse_complement(id = assign_result.qid, description = "reverse")
			assign_position_dict[assign_result.qid] = [origin_seq.seq, assign_result.qstart, assign_result.qend, assign_result.sstart]
		if assign_result.assign_type == "J":
			try:
				assign_position_dict[assign_result.qid].extend([assign_result.qend, assign_result.qstart])
			except ValueError:
				pass

	for key, value in assign_position_dict.items():
		query_seq, start, V_gene_region_end, trans_start = value[0], int(value[1]), int(value[2]), int(value[3])
		try:
			Variable_region_end = int(value[4])
			J_start = int(value[5])
			result[7] = query_seq[J_start : Variable_region_end]
			J_gene_region   = SeqRecord_gernerator(key, str(result[7]), 'J_gene_region')
			SeqIO.write(J_gene_region, outfile3, "fasta")
		except IndexError:
			continue
		if (trans_start-1) % 3 == 0:
			result[3] = query_seq[start -1 :Variable_region_end]
			result[4] = query_seq[start -1 :V_gene_region_end]
			result[5] = query_seq[ : start -1]
			result[6] = query_seq[Variable_region_end :]
 		else:
			result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):Variable_region_end]
			result[4] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):V_gene_region_end]
			result[5] = query_seq[ : (start - 1 + 2 - ((trans_start-1)%3 -1))]
			result[6] = query_seq[Variable_region_end :]
		try:
			CDR3_region_part = query_seq[int(CDR3_start_dict[key]) - 1 : ]
			#CDR3_region_part_protein = CDR3_region_part.translate()
			has_fwgxg = re.search('((TGG)|(TTT)|(TTC))GG[AGCT][AGCT][AGCT][AGCT]GG[AGCT]',str(CDR3_region_part))
			try:
				CDR3_region = CDR3_region_part[ : has_fwgxg.end()]
				CDR3_region_out = SeqRecord_gernerator(key, str(CDR3_region), 'CDR3_region')
			except:
				CDR3_region = query_seq[int(CDR3_start_dict[key]) - 1 : J_start]
				CDR3_region_out = SeqRecord_gernerator(key, str(CDR3_region), 'CDR3_region_no_wgxg')
			SeqIO.write(CDR3_region_out, CDR3_outfile, "fasta")
		except KeyError:
			pass
		Variable_region = SeqRecord_gernerator(key, str(result[3]), 'Variable_region')
		V_gene_region   = SeqRecord_gernerator(key, str(result[4]), 'V_gene_region')
		SeqIO.write(Variable_region, outfile, "fasta")
		SeqIO.write(V_gene_region, outfile2, "fasta")
		if result[5] != '':
			Leader_region = SeqRecord_gernerator(key, str(result[5]), 'Leader_region')
			SeqIO.write(Leader_region, leader_region_outfile, "fasta")
		if result[6] != '':
			C_region = SeqRecord_gernerator(key, str(result[6]), 'C_region')
			SeqIO.write(C_region, C_region_outfile, "fasta")
	CDR3_outfile.close()
	outfile.close()
	outfile2.close()
	outfile3.close()
	C_region_outfile.close()
	leader_region_outfile.close()
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
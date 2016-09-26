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
def trim_Variable_region(prj_tree, prj_name, IGBLAST_assignment_file):
	outfile = open("%s/%s_trimmed_Variable_region.fasta"%(prj_tree.reads, prj_name),"w")
	outfile2 = open("%s/%s_trimmed_V_gene_region.fasta"%(prj_tree.reads, prj_name),"w")
	record_dict = SeqIO.index('%s/%s.assembled_trimed.fasta'%(prj_tree.origin, prj_name),  "fasta")
	reader = csv.reader(open(IGBLAST_assignment_file,"rU"), delimiter = "\t")
	result = ['Query_ID','query_seq','query_length','Variable_region', 'V_gene_region']
	ID = ''
	count = 0
	count1 = 0
	start, Variable_region_end, V_gene_region = int(), int(), int()
	for line in reader:
		if line !=[]:
			query_ID = line[1].replace(' ','')
			nr_query_ID = query_ID.replace('reversed|','')
			if query_ID != ID and line[0].replace(' ','') == "V":
				count1 += 1
				if (count1) % 10000 == 0:
					print "%d reads processed" %(count1)
				#print nr_query_ID
				query_seq = str(record_dict.get(nr_query_ID).seq)
				if 'reversed' in query_ID:
					count += 1
					#if (count) % 1000 == 0:
						#print "%d reversed reads processed" %(count)
					a = record_dict.get(nr_query_ID)
					rc = a.reverse_complement(id = query_ID )
					query_seq = str(rc.seq)
				#writer.writerow(result)
				if result[0] != 'Query_ID':
					Variable_region = SeqRecord_gernerator(result[0], result[3], 'Variable_region')
					V_gene_region = SeqRecord_gernerator(result[0], result[4], 'V_gene_region')
					SeqIO.write(Variable_region, outfile, "fasta")
					SeqIO.write(V_gene_region, outfile2, "fasta")
				ID = query_ID
				start = int(line[8])
				#print line[10],line
				trans_start = int(line[10])
				V_gene_region_end = int(line[9])
			elif line[0].replace(' ','') == "J":
				Variable_region_end = int(line[9])
			result[0] = nr_query_ID
			result[1] = query_seq
			result[2] = len(query_seq)
			if start-trans_start < 0:
				
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :Variable_region_end]
					result[4] = query_seq[start -1 :V_gene_region_end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):Variable_region_end]
					result[4] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):V_gene_region_end]
			else:
				if (trans_start-1) % 3 == 0:
					result[3] = query_seq[start -1 :Variable_region_end]
					result[4] = query_seq[start -1 :V_gene_region_end]
				else:
					#Wrong! and beside correceted! result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
					result[3] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):Variable_region_end]
					result[4] = query_seq[(start - 1 + 2 - ((trans_start-1)%3 -1)):V_gene_region_end]
	#writer.writerow(result)
	Variable_region = SeqRecord_gernerator(result[0], result[3], 'Variable_region')
	V_gene_region = SeqRecord_gernerator(result[0], result[4], 'V_gene_region')
	SeqIO.write(Variable_region, outfile, "fasta")
	SeqIO.write(V_gene_region, outfile2, "fasta")

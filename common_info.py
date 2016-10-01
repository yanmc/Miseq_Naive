#!/usr/bin/env python
# encoding: utf-8
"""
common_info.py

Created by Zhenhai Zhang on 2011-04-06.
Copyright (c) 2011 Columbia University and Vaccine Research Center, National Institutes of Health, USA. All rights reserved.
"""

ORG_FOLDER = "origin"
FILTERED_FOLDER="filtered"
MAPPING_FOLDER="mapping" 
ANALYSIS_FOLDER="analysis"
LOG_FOLDER="logs"
CLUSTAL_FOLDER="clustal"
PHYLO_FOLDER="phylo"
DOC_FOLDER="docs"
TMP_FOLDER="tmp"

SPLIT_FOLDER="split" 
GERM_FOLDER="germ"
NAT_FOLDER="native"	
SELF_FOLDER="reads"
PBS_FOLDER="pbs"
JOB_FOLDER="jobs"

ANALYSIS_DATA_FOLDER="data"
ANALYSIS_FIGURE_FOLDER="figure"

CLUSTAL_FASTA_FOLDER="clustal_fasta"
CLUSTAL_PBS_FOLDER="clustal_pbs"
CLUSTAL_JOB_FOLDER="clustal_job"
CLUSTAL_DATA_FOLDER ="clustal_data"

IGBLAST_DATA_FOLDER = "Igblast_data"
IGBLAST_DATABASE_FOLDER = "Igblast_database"
IGBLAST_TMP_FOLDER = "Igblast_tmp"


nucleotides = ["A", "C", "G", "T"]

#
# ===START=== amino acid
#
dict_codon2aa = {

	"ATT" : ("I", "Isoleucine"),
	"ATA" : ("I", "Isoleucine"),
	"ATC" : ("I", "Isoleucine"),
	
	"CTT" : ("L", "Leucine"),
	"CTC" : ("L", "Leucine"),
	"CTA" : ("L", "Leucine"),
	"CTG" : ("L", "Leucine"),
	"TTA" : ("L", "Leucine"),
	"TTG" : ("L", "Leucine"),
	
	"GTT" : ("V", "Valine"),
	"GTC" : ("V", "Valine"),
	"GTA" : ("V", "Valine"),
	"GTG" : ("V", "Valine"),
	
	"TTT" : ("F", "Phenylalanine"),
	"TTC" : ("F", "Phenylalanine"),
	
	"ATG" : ("M", "Methionine"),
	
	"TGT" : ("C", "Cysteine"),
	"TGC" : ("C", "Cysteine"),
	
	"GCT" : ("A", "Alanine"),
	"GCC" : ("A", "Alanine"),
	"GCA" : ("A", "Alanine"),
	"GCG" : ("A", "Alanine"),
	
	"GGT" : ("G", "Glycine"),
	"GGC" : ("G", "Glycine"),
	"GGA" : ("G", "Glycine"),
	"GGG" : ("G", "Glycine"),
	
	"CCT" : ("P", "Proline"),
	"CCC" : ("P", "Proline"),
	"CCA" : ("P", "Proline"),
	"CCG" : ("P", "Proline"),
	
	"ACT" : ("T", "Threonine"),
	"ACC" : ("T", "Threonine"),
	"ACA" : ("T", "Threonine"),
	"ACG" : ("T", "Threonine"),
	
	"TCT" : ("S", "Serine"),
	"TCC" : ("S", "Serine"),
	"TCA" : ("S", "Serine"),
	"TCG" : ("S", "Serine"),
	"AGT" : ("S", "Serine"),
	"AGC" : ("S", "Serine"),
	
	"TAT" : ("Y", "Tyrosine"),
	"TAC" : ("Y", "Tyrosine"),
	
	"TGG" : ("W", "Tryptophan"),
	
	"CAA" : ("Q", "Glutamine"),
	"CAG" : ("Q", "Glutamine"),
	
	"AAT" : ("N", "Asparagine"), 
	"AAC" : ("N", "Asparagine"),
	
	"CAT" : ("H", "Histidine"), 
	"CAC" : ("H", "Histidine"),
	
	"GAA" : ("E", "Glutamic acid"),
	"GAG" : ("E", "Glutamic acid"),
	
	"GAT" : ("D", "Aspartic acid"),
	"GAC" : ("D", "Aspartic acid"),
	
	"AAA" : ("K", "Lysine"),
	"AAG" : ("K", "Lysine"),
	
	"CGT" : ("R", "Arginine"),
	"CGC" : ("R", "Arginine"),
	"CGA" : ("R", "Arginine"),
	"CGG" : ("R", "Arginine"),
	"AGA" : ("R", "Arginine"),
	"AGG" : ("R", "Arginine"),
	
	"TAA" : ("STOP", "Stop codons"),
	"TAG" : ("STOP", "Stop codons"),
	"TGA" : ("STOP", "Stop codons")

	}
	
start_codons = {"ATG"}
stop_codons = {"TAA", "TAG", "TGA"}

dict_aa2codon = {

	"I" : ("ATT", "ATC", "ATA"),
	"L" : ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
	"V" : ("GTT", "GTC", "GTA", "GTG"),
	"F" : ("TTT", "TTC"),
	"M" : ("ATG"),
	"C" : ("TGT", "TGC"),
	"A" : ("GCT", "GCC", "GCA", "GCG"),
	"G" : ("GGT", "GGC", "GGA", "GGG"),
	"P" : ("CCT", "CCC", "CCA", "CCG"),
	"T" : ("ACT", "ACC", "ACG", "ACA"),
	"S" : ("TCT", "TCC", "TCG", "TCA", "AGT", "AGC"),
	"Y" : ("TAT", "TAC"),
	"W" : ("TGG"),
	"Q" : ("CAA", "CAG"),
	"N" : ("AAT", "AAC"),
	"H" : ("CAT", "CAC"),
	"E" : ("GAA", "GAG"),
	"D" : ("GAT", "GAC"),
	"K" : ("AAA", "AAG"),
	"R" : ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
	"STOP" : ("TAA", "TAG", "TGA")

	}
dict_aa2name = {

	"I" : "Isoleucine",
	"L" : "Leucine",
	"V" : "Valine",
	"F" : "Phenylalanine",
	"M" : "Methionine",
	"C" : "Cysteine",
	"A" : "Alanine",
	"G" : "Glycine",
	"P" : "Proline",
	"T" : "Threonine",
	"S" : "Serine",
	"Y" : "Tyrosine",
	"W" : "Tryptophan",
	"Q" : "Glutamine",
	"N" : "Asparagine",
	"H" : "Histidine",
	"E" : "Glutamic acid",
	"D" : "Aspartic acid",
	"K" : "Lysine",
	"R" : "Arginine",
	"STOP" : "Stop condons"
}

#
# ===END=== amino acid
#


#
# ===START=== amino acid
#
IUB_CODE = {
	"A" : ["A"], 
	"C" : ["C"], 
	"G" : ["G"], 
	"T" : ["T"], 
	"R" : ["A", "G"], 
	"Y" : ["C", "T"], 
	"K" : ["G", "T"], 
	"M" : ["A", "C"], 
	"S" : ["G", "C"], 
	"W" : ["A", "T"], 
	"B" : ["C", "G", "T"], 
	"D" : ["A", "G", "T"], 
	"H" : ["A", "C", "T"], 
	"V" : ["A", "C", "G"], 
	"N" : ["A", "C", "G", "T"]
}
	
IUB_complement = {
	"A" : "T",
	"C" : "G", 
	"G" : "C", 
	"T" : "A",
	"R" : "Y", 
	"Y" : "R",
	"K" : "M", 
	"M" : "K",
	"S" : "W",
	"W" : "S",
	"B" : "V",
	"V" : "B",
	"D" : "H",
	"H" : "D",
	"N" : "N"
}
#
# ===END=== amino acid
#


#
# ===START=== HUMAN GERMLINE
#
HUMAN_GERMLINE = {
"HUMANIGHV" : ['IGHV4-4', 'IGHV3-19', 'IGHV3-15', 'IGHV1-3', 'IGHV1-18', 'IGHV3-11', 'IGHV3-13', 'IGHV3-33-2', 'IGHV3-7', 'IGHV3-9', 'IGHV3-30-3', 'IGHV3-30-5', 'IGHV1-24', 'IGHV3-16', 'IGHV1-NL1', 'IGHV4-30-4', 'IGHV2-70', 'IGHV4-30-2', 'IGHV2-70D', 'IGHV3-62', 'IGHV3-63', 'IGHV3-64', 'IGHV3-66', 'IGHV7-4-1', 'IGHV5-78', 'IGHV7-81', 'IGHV3/OR15-7', 'IGHV3/OR16-9', 'IGHV3/OR16-8', 'IGHV3-73', 'IGHV3-72', 'IGHV3-71', 'IGHV3/OR16-6', 'IGHV3-74', 'IGHV1-8', 'IGHV1-38-4', 'IGHV1-2', 'IGHV4-39', 'IGHV4-34', 'IGHV4-38-2', 'IGHV4-31', 'IGHV4/OR15-8', 'IGHV3-48', 'IGHV3-49', 'IGHV3-47', 'IGHV3-43', 'IGHV3-30-22', 'IGHV1-69-2', 'IGHV3-30-2', 'IGHV3-NL1', 'IGHV4-28', 'IGHV1/OR21-1', 'IGHV6-1', 'IGHV3-30-33', 'IGHV3-53', 'IGHV3-52', 'IGHV3-54', 'IGHV2-5', 'IGHV4-55', 'IGHV3-69-1', 'IGHV1-45', 'IGHV1-46', 'IGHV4-59', 'IGHV3-25', 'IGHV3-20', 'IGHV3-30-42', 'IGHV3-22', 'IGHV3-23', 'IGHV7-34-1', 'IGHV3-29', 'IGHV7-40', 'IGHV1-58', 'IGHV3-33', 'IGHV5-51', 'IGHV3/OR16-15', 'IGHV3/OR16-14', 'IGHV3/OR16-16', 'IGHV3/OR16-10', 'IGHV3/OR16-13', 'IGHV3/OR16-12', 'IGHV3-23D', 'IGHV3-35', 'IGHV1-69D', 'IGHV3-32', 'IGHV3-30-52', 'IGHV3-30', 'IGHV2-26', 'IGHV3-38', 'IGHV5-10-1', 'IGHV3-64D', 'IGHV1-68', 'IGHV1-69', 'IGHV3-38-3', 'IGHV1/OR15-1', 'IGHV1/OR15-2', 'IGHV1/OR15-3', 'IGHV1/OR15-4', 'IGHV1/OR15-5', 'IGHV1/OR15-9', 'IGHV2-10', 'IGHV4-61', 'IGHV2/OR16-5', 'IGHV3-43D', 'IGHV3-21'],
"HUMANIGHD" : ['IGHD4-23', 'IGHD6-6', 'IGHD1-20', 'IGHD6-19', 'IGHD1-26', 'IGHD6-13', 'IGHD3-22', 'IGHD1/OR15-1a', 'IGHD1/OR15-1b', 'IGHD2-15', 'IGHD2-2', 'IGHD1-7', 'IGHD1-1', 'IGHD2-8', 'IGHD5-24', 'IGHD5/OR15-5a', 'IGHD5/OR15-5b', 'IGHD4-17', 'IGHD1-14', 'IGHD4-11', 'IGHD4-4', 'IGHD3-9', 'IGHD5-5', 'IGHD3/OR15-3a', 'IGHD3/OR15-3b', 'IGHD3-3', 'IGHD2-21', 'IGHD6-25', 'IGHD3-16', 'IGHD3-10', 'IGHD7-27', 'IGHD4/OR15-4b', 'IGHD5-18', 'IGHD4/OR15-4a', 'IGHD5-12', 'IGHD2/OR15-2b', 'IGHD2/OR15-2a'],
"HUMANIGHJ" : ['IGHJ6', 'IGHJ5', 'IGHJ4', 'IGHJ3', 'IGHJ2', 'IGHJ1'],
"HUMANIGKV" : ['IGKV4-1', 'IGKV1/OR2-11', 'IGKV3/OR2-268', 'IGKV1-12', 'IGKV1-13', 'IGKV1/OR2-9', 'IGKV1-16', 'IGKV1-17', 'IGKV2D-30', 'IGKV2-4', 'IGKV1D-12', 'IGKV1/OR2-108', 'IGKV1-37', 'IGKV1/OR1-1', 'IGKV2-24', 'IGKV1-33', 'IGKV1D-17', 'IGKV1D-16', 'IGKV1D-33', 'IGKV2-29', 'IGKV1/OR2-2', 'IGKV1-39', 'IGKV1/OR2-0', 'IGKV1/OR2-1', 'IGKV2-40', 'IGKV5-2', 'IGKV6-21', 'IGKV1D-13', 'IGKV1/ORY-1', 'IGKV1-NL1', 'IGKV2-28', 'IGKV1/OR2-3', 'IGKV3D-20', 'IGKV3-11', 'IGKV1D-37', 'IGKV2D-24', 'IGKV1/OR10-1', 'IGKV6D-41', 'IGKV2D-29', 'IGKV2D-28', 'IGKV1/OR22-5', 'IGKV3-7', 'IGKV2D-18', 'IGKV2D-26', 'IGKV1/OR9-1', 'IGKV1D-42', 'IGKV1/OR9-2', 'IGKV1-27', 'IGKV3-20', 'IGKV1/OR2-118', 'IGKV2/OR2-7D', 'IGKV2-30', 'IGKV2-18', 'IGKV3D-7', 'IGKV1D-43', 'IGKV1D-39', 'IGKV1D-8', 'IGKV1-5', 'IGKV1-6', 'IGKV1-8', 'IGKV1-9', 'IGKV3-15', 'IGKV2/OR22-4', 'IGKV6D-21', 'IGKV1/OR15-118', 'IGKV3D-11', 'IGKV7-3', 'IGKV1/OR-3', 'IGKV1/OR-2', 'IGKV3D-15', 'IGKV1/OR-4', 'IGKV2D-40'],
"HUMANIGKJ" : ['IGKJ1', 'IGKJ2', 'IGKJ3', 'IGKJ4', 'IGKJ5'],
"HUMANIGLV" : ['IGLV2-8', 'IGLV4-3', 'IGLV5-45', 'IGLV1-62', 'IGLV5-48', 'IGLV2-5', 'IGLV2-11', 'IGLV6-57', 'IGLV2-14', 'IGLV2-34', 'IGLV9-49', 'IGLV2-18', 'IGLV1-36', 'IGLV2-33', 'IGLV7-46', 'IGLV3-10', 'IGLV3-13', 'IGLV3-12', 'IGLV7-43', 'IGLV3-16', 'IGLV3-19', 'IGLV3-32', 'IGLV1-51', 'IGLV1-50', 'IGLV3-27', 'IGLV11-55', 'IGLV8-61', 'IGLV5-52', 'IGLV5-37', 'IGLV5-39', 'IGLV8/OR8-1', 'IGLV10-54', 'IGLV2-23', 'IGLV2-NL1', 'IGLV1-40', 'IGLV1-41', 'IGLV3-22', 'IGLV1-44', 'IGLV3-25', 'IGLV1-47', 'IGLV4-60', 'IGLV3-1', 'IGLV3-21', 'IGLV4-69', 'IGLV3-9', 'IGLV3-31'],
"HUMANIGLJ" : ['IGLJ3', 'IGLJ2', 'IGLJ1', 'IGLJ7', 'IGLJ6', 'IGLJ5', 'IGLJ4']
}
#
# ===END=== HUAMN GERMLINE
#

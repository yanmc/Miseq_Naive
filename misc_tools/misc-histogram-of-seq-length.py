#!/usr/bin/env python
# encoding: utf-8
"""
8.0-histogram-of-seq-length.py -i infile -org orgnism

Created by Mingchen on 2014-12-16.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

orgnism:  human,mouse or rabbit
infile: trimmed_seq: the file contain trimmed sequence in 6.41-get-trimmed-reads-fasta

"""
from mytools import *
from Bio import SeqIO
import pylab

import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import pandas as pd
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
#from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
#from misc_prepare_pbs import *
#from misc_get_trimmed_region import *
from collections import Counter
try:
    import cPickle as pickle
except ImportError:
    import pickle
import statsmodels.api as sm
def main():
	#infile = "%s/%s.assembled.fastq"%(prj_tree.origin, prj_name)
	infile = "%s/%s_H_leader_region.fasta"%(prj_tree.reads, prj_name)
	name, suffix = os.path.splitext(infile)
	sizes = [len(rec) for rec in SeqIO.parse(infile, "fasta")]
	print len(sizes), min(sizes), max(sizes)
	#print sizes
	pylab.hist(sizes, bins=20)
	pylab.title("%i orchid sequences\nLengths %i to %i" % (len(sizes),min(sizes),max(sizes)))
	pylab.xlabel("Sequence length (bp)")
	pylab.ylabel("count")
	#pylab.show()
	pylab.savefig('%s/%s_Leader_histogram_length.png'%(prj_tree.figure, prj_name)) #to save the figure to a file (e.g. as a PNG or PDF).
if __name__ == '__main__':
	
	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()

	main()
	end = time.time()
	print prj_name, end-start
	print "Finished"
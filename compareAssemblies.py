#!/usr/bin/env python3

import mappy as mp
import pandas as pd
import csv
import argparse
import re
import numpy as np
from Bio import SeqIO,SeqUtils
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

# Notes: gzip not supported with biopython SeqIO

def parseArgs() :
	parser = argparse.ArgumentParser(description='compare assemblies with minimap2')
	parser.add_argument('-r', help='input assembly 1 (ref) fasta', dest='asm1Filename',required=True)
	parser.add_argument('-q', help='input assembly 2 (query) fasta', dest='asm2Filename',required=True)
	parser.add_argument('-m', help='minimum gap size, default=25', dest='minGapSize',required=False, default=25, type=int)
	parser.add_argument('-n', help='minimum sequence length in assembly 2 (query), default=20Mb', dest='minQueryLen',required=False, default=20000000, type=int)
	parser.add_argument('-o', help='output prefix, default=output', dest='outputPrefix', required=False, default='output')
	parser.add_argument('-d', help='turn off saving minimap2 results', dest='saveMinimap', required=False, action='store_false', default=True)
	arguments = parser.parse_args()
	print("arguments:")
	print("-r, input assembly 1 (ref) fasta: %s" % arguments.asm1Filename)
	print("-q, input assembly 2 (query) fasta: %s" % arguments.asm2Filename)
	print("-m, minimum gap size: %i" % arguments.minGapSize)
	print("-n, minimum sequence length in assembly 2: %i" % arguments.minQueryLen)
	print("-o, output prefix: %s" % arguments.outputPrefix)
	print("(-d), save minimap2 results?: %s" % arguments.saveMinimap)
	print("-- \n")
		
	return arguments

def runIndex(refFastaFilname) :
	print("building index for assembly 1 (ref): %s\n" % refFastaFilname)
	a = mp.Aligner(refFastaFilname, preset="asm5")
	return a


def getTopHitByAlignmentLength(hitsList, queryID) :
	"function getTopHitByAlignmentLength"
	hitsDf = pd.DataFrame([i.split('\t') for i in hitsList], columns=['q','q_len','q_st','q_en','strand','ctg','ctg_len','r_st','r_en','mlen','blen','mapq','f1','f2','cigar'])
	hitsDf['blen'] = hitsDf['blen'].astype(int)
	hitsDf_subset = hitsDf[hitsDf['q']==queryID]
	grouped = hitsDf_subset.groupby('ctg')
	# get ctg with most total alignments
	if grouped :
		output = {'top_aln_id': grouped['blen'].sum().idxmax(), \
              'top_aln_blen': grouped['blen'].sum().max()
              }
	else :
		output = {'top_aln_id': 'NONE', \
              'top_aln_blen': 'NONE'
              }
	return output

def runMapper(referenceIndex, asm2Filename, minQueryLen, saveMinimap) :
	print("running minimap2 and finding top hit per query sequence\n")             
	scaffoldMapList0 = []
	hitsListAll0 = []
	for name, seq, qual in mp.fastx_read(asm2Filename):
		print("... query: %s" % name)
		if len(seq) < minQueryLen:
			print("...... Skipping, query too short (seq len of %i is less than minimum: %i)\n" % (len(seq), minQueryLen))
			continue
		hits=[]                 
		for hit in referenceIndex.map(seq):
			if hit.is_primary:
				hits.append(name+"\t"+str(len(seq))+"\t"+str(hit))
		if hits:
			if saveMinimap:
					hitsListAll0.append(name+"\t"+str(len(seq))+"\t"+str(hit))
			topAln = getTopHitByAlignmentLength(hits, name)
			print("Top hit: %s\n" % topAln['top_aln_id'])        
			scaffoldMapList0.append({'queryID': name, 
								 'refID' : topAln['top_aln_id'], 
								 'alignLen' : topAln['top_aln_blen'],
								 })
	return (scaffoldMapList0, hitsListAll0)

def writeOutput(outputPrefix, saveMinimap, scaffoldMapListOut, hitsListAllOut) :
	print("writing results")
	colnames = list(scaffoldMapListOut[0].keys())
	colnames.sort()
	with open(outputPrefix + "_top_hit_results.csv", 'w') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames= colnames) 
		writer.writeheader()
		for data in scaffoldMapListOut:
			writer.writerow(data)
	
	if saveMinimap:
		hitsDf = pd.DataFrame([i.split('\t') for i in hitsListAllOut], columns=['q','q_len','q_st','q_en','strand','ctg','ctg_len','r_st','r_en','mlen','blen','mapq','f1','f2','cigar'])
		hitsDf.to_csv(outputPrefix + "_minimap2.txt", index=False, sep="\t")

def getMetrics(scaffoldMapList, asm1, asm2, min_gap_size, outputPrefix):
	ref_dict = SeqIO.index(asm1, "fasta")
	query_dict = SeqIO.index(asm2, "fasta")
	ref_lengths = []
	query_lengths = []
	# loop over scaffold pairs, get info on lengths, gaps
	for i in scaffoldMapList:
		# reference asm gaps
		gaps_ref = re.findall('([Nn]+)', str(ref_dict[i['refID']].seq))
		gaps_ref_filtered = [x for x in gaps_ref if len(x) >= min_gap_size]
		i['ref_numGaps'] = len(gaps_ref_filtered) # append new data to dictionary scaffoldMapList
		i['ref_totalGapsLength'] = len("".join(gaps_ref_filtered))  # append new data to dictionary scaffoldMapList
		i['ref_meanGapsLength'] = np.mean([len(x) for x in gaps_ref_filtered]).round(decimals=2)
		i['ref_length'] = len(ref_dict[i['refID']].seq)
		ref_lengths.append(i['ref_length'])
	
		# query asm gaps
		gaps_query = re.findall('([Nn]+)', str(query_dict[i['queryID']].seq))
		gaps_query_filtered = [x for x in gaps_query if len(x) >= min_gap_size]
		i['query_numGaps'] = len(gaps_query_filtered)
		i['query_totalGapsLength'] = len("".join(gaps_query_filtered))
		i['query_meanGapsLength'] = np.mean([len(x) for x in gaps_query_filtered]).round(decimals=2)
		i['query_length'] = len(query_dict[i['queryID']].seq)
		query_lengths.append(i['query_length'])
		
		# GC content
		i['ref_GC'] = SeqUtils.GC(ref_dict[i['refID']].seq)
		i['query_GC'] = SeqUtils.GC(query_dict[i['queryID']].seq)
	
	# scaffold length differences
	lengthDiffs = [x['ref_length']-x['query_length'] for x in scaffoldMapList]
	gapLengthDiffs = [x['ref_totalGapsLength']-x['query_totalGapsLength'] for x in scaffoldMapList]
	gapCountDiffs = [x['ref_numGaps']-x['query_numGaps'] for x in scaffoldMapList]

	figSizes = plt.figure()
	sns.distplot(lengthDiffs, bins=20)
	figSizes.savefig(outputPrefix + "_histogram-scaffold_diffs.pdf", bbox_inches='tight')
	plt.close()
	
	# n50 - calc by slicing list of scaffolds lengths, using the index at the n50 entry
	ref_lengths.sort(reverse=True)
	query_lengths.sort(reverse=True)
	ref_n50 = ref_lengths[len([x for x in np.cumsum(ref_lengths) if x <= np.sum(ref_lengths)*.5])-1]
	query_n50 = query_lengths[len([x for x in np.cumsum(query_lengths) if x <= np.sum(query_lengths)*.5])-1]
	
	# write summary metrics
	summaryOut = open(outputPrefix + "_summary.txt", "w")
	summaryOut.write("Comparison between between scaffold pairs, 'Ref-value' minus 'Query-value' \n")
	summaryOut.write("Average difference in scaffold lengths: %s bp\n" % np.mean(lengthDiffs).round(decimals=2))
	summaryOut.write("Average difference in gap content: %s bp\n" % np.mean(gapLengthDiffs).round(decimals=2))
	summaryOut.write("Average difference in number of gaps: %s\n\n" % np.mean(gapCountDiffs).round(decimals=2))

	summaryOut.write("Total difference in scaffold lengths: %s bp\n" % np.sum(lengthDiffs))
	summaryOut.write("Total difference in gap content: %s bp\n" % np.sum(gapLengthDiffs))
	summaryOut.write("Total difference in number of gaps: %s\n\n" % np.sum(gapCountDiffs))
	
	summaryOut.write("Ref scaffolds N50: %s\n" % ref_n50)
	summaryOut.write("Query scaffolds N50: %s\n" % query_n50)

	summaryOut.close()

	return scaffoldMapList


def run() :
	"run script"
	# parse arguments
	args = parseArgs()
	
	# generate index for reference fasta
	refIndex = runIndex(refFastaFilname = args.asm1Filename)
	
	# run minimap2 - loop over query sequences
	(scaffoldMapList, hitsListAll) = runMapper(referenceIndex = refIndex, asm2Filename = args.asm2Filename, minQueryLen = args.minQueryLen, saveMinimap = args.saveMinimap)

	# run gap finding
	scaffoldMapList = getMetrics(scaffoldMapList=scaffoldMapList, asm1=args.asm1Filename, asm2=args.asm2Filename, min_gap_size = args.minGapSize, outputPrefix = args.outputPrefix)

	# write output
	writeOutput(outputPrefix = args.outputPrefix, saveMinimap = args.saveMinimap, scaffoldMapListOut = scaffoldMapList, hitsListAllOut = hitsListAll)
	

	print("DONE")

if __name__ == '__main__':
    run()


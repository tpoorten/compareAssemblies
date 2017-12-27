#!/usr/bin/env python3

import mappy as mp
import pandas as pd
import csv
import argparse

def parseArgs() :
	parser = argparse.ArgumentParser(description='compare assemblies with minimap2')
	parser.add_argument('-a', help='input assembly 1 (ref) fasta', dest='asm1Filename',required=True)
	parser.add_argument('-b', help='input assembly 2 (query) fasta', dest='asm2Filename',required=True)
	parser.add_argument('-m', help='minimum sequence length in assembly 2, default=20Mb', dest='minQueryLen',required=False, default=20000000, type=int)
	parser.add_argument('-o', help='output prefix, default=output', dest='outputPrefix', required=False, default='output')
	parser.add_argument('-d', help='turn off saving minimap2 results', dest='saveMinimap', required=False, action='store_false', default=True)
	arguments = parser.parse_args()
	print("arguments:")
	print("input assembly 1 (ref) fasta: %s" % arguments.asm1Filename)
	print("input assembly 2 (query) fasta: %s" % arguments.asm2Filename)
	print("minimum sequence length in assembly 2: %i" % arguments.minQueryLen)
	print("output prefix: %s" % arguments.outputPrefix)
	print("turn off saving minimap2 results: %s" % arguments.saveMinimap)
	print("-- \n")
		
	return arguments

def runIndex(refFastaFilname) :
	print("building index for assembly 1 (ref): %s\n" % refFastaFilname)
	a = mp.Aligner(refFastaFilname, preset="asm5")
	return a


def getTopHitByAlignmentLength(hitsList, queryID) :
	"function getTopHitByAlignmentLength"
	hitsDf = pd.DataFrame([i.split('\t') for i in hitsList], columns=['q','q_len','q_st','q_en','strand','ctg','ctg_len','r_st','r_en','mlen','blen','mapq','f1','f2','cigar'])
	hitsDf['ctg_len'] = hitsDf['ctg_len'].astype(int)
	hitsDf_subset = hitsDf[hitsDf['q']==queryID]
	grouped = hitsDf_subset.groupby('ctg')
	# get ctg with most total alignments
	output = {'top_aln_id': grouped['ctg_len'].sum().idxmax(), \
              'top_aln_len':grouped['ctg_len'].sum().max()}
	return output

def runMapper(referenceIndex, asm2Filename, minQueryLen, saveMinimap) :
	print("running minimap2 and finding top hit per query sequence\n")             
	scaffoldMapList0 = []
	hitsListAll0 = []
	for name, seq, qual in mp.fastx_read(asm2Filename):
		print("... query: %s" % name)
		if len(seq) < minQueryLen:
			print("...... Skipping, query too short (seq len of %i is less than minimum: %i)\n" % (len(seq), args.minQueryLen))
			continue
		hits=[]                 
		for hit in referenceIndex.map(seq):
			if hit.is_primary:
				hits.append(name+"\t"+str(len(seq))+"\t"+str(hit))
				if saveMinimap:
					hitsListAll0.append(name+"\t"+str(len(seq))+"\t"+str(hit))
		topAln = getTopHitByAlignmentLength(hits, name)
		print("Top hit: %s\n" % topAln['top_aln_id'])        
		scaffoldMapList0.append({'queryID': name, 'refID' : topAln['top_aln_id'], 'alignLen' : topAln['top_aln_len']})
	return (scaffoldMapList0, hitsListAll0)

def writeOutput(outputPrefix, saveMinimap, scaffoldMapListOut, hitsListAllOut) :
	print("writing results")
	with open(outputPrefix + "_top_hit_results.csv", 'w') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=['queryID','refID','alignLen'])
		writer.writeheader()
		for data in scaffoldMapListOut:
			writer.writerow(data)
	
	if saveMinimap:
		hitsDf = pd.DataFrame([i.split('\t') for i in hitsListAllOut], columns=['q','q_len','q_st','q_en','strand','ctg','ctg_len','r_st','r_en','mlen','blen','mapq','f1','f2','cigar'])
		hitsDf.to_csv(outputPrefix + "_minimap2.txt", index=False, sep="\t")


def run() :
	"run script"
	# parse arguments
	args = parseArgs()
	
	# generate index for reference fasta
	refIndex = runIndex(refFastaFilname = args.asm1Filename)
	
	# run minimap2 - loop over query sequences
	(scaffoldMapList, hitsListAll) = runMapper(referenceIndex = refIndex, asm2Filename = args.asm2Filename, minQueryLen = args.minQueryLen, saveMinimap = args.saveMinimap)

	# write output
	writeOutput(outputPrefix = args.outputPrefix, saveMinimap = args.saveMinimap, scaffoldMapListOut = scaffoldMapList, hitsListAllOut = hitsListAll)

	print("DONE")

if __name__ == '__main__':
    run()


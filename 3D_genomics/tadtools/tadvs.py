## calculate intra TAD vs inter TAD
## consolidation score

import pandas as pd
import argparse
import sys
import statistics

parser = argparse.ArgumentParser()

parser.add_argument('-m', type=str, help='dense matrix to analyze')
parser.add_argument('-t', type=str, help='tads file as reference')

args = parser.parse_args()

def getTAD(tadfile):
	tads = pd.read_csv(tadfile,sep='\t',header=0,usecols=[0,1,2],names=['chr','start','end'],dtype={'chr':'str','start':'int','end':'int'})
	return tads

def getSize(tads):
	sizes = dict.fromkeys(['chr1','chr2','chr3','chr4','chr5'])
	for chrom in sizes:
		sizes[chrom] = tads.loc[tads['chr']==chrom].iloc[-1]['end']
	return sizes

def getMat(matfile):
	mat = pd.read_csv(matfile,sep='\t',header=0,names=['chr','first','sec','value'],dtype={'chr':'str','first':'int','sec':'int','value':'float'})
	## remove short-distance region (<100kb)
	#mat = mat.loc[mat['sec']-mat['first']>100000]
	return mat

def getValue(tad,mat,sizes):
	length = sizes[tad[0]]
	## select chr
	matChr = mat.loc[mat['chr']==tad[0]]
	## select region
	matRe = mat.loc[(mat['sec']>mat['first']) & (mat['sec']<mat['first']+tad[2]-tad[1]) & (mat['sec']>2*tad[1]-mat['first']) & (mat['sec']<2*tad[2]-mat['first'])]
	matIntra = matRe.loc[(matRe['first']>=tad[1]) & (matRe['sec']<=tad[2])]
	matInter = matRe.loc[(matChr['first']<tad[1]) | (matRe['sec']>tad[2])]
	if (matIntra.shape[0]>0) & (matInter.shape[0]>0):
		value = statistics.mean(matIntra['value'])/statistics.mean(matInter['value'])
		return value
	else:
		return 'not'

if __name__ == '__main__':

	tads = getTAD(args.t)
	mat = getMat(args.m)
	sizes = getSize(tads)

	sys.stdout.write("[%s]" % (" " * tads.shape[0]))
	sys.stdout.flush()
	sys.stdout.write("\b" * (tads.shape[0]+1))

	with open(args.m+'_all.txt','a') as f:
		for row in tads.itertuples(index=False):
			ci = getValue(row,mat,sizes)
			if ci != 'not':
				f.write('%s\n' %ci)
			sys.stdout.write('-')
			sys.stdout.flush()

	sys.stdout.write(']\n')








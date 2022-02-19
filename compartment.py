import os,sys
import argparse

parser = argparse.ArgumentParser(description="PCA analysis for C-data \
using HOMER")
parser.add_argument("-t", "--makeTagDirectory", action="store_true",
    default=False)
parser.add_argument("-r", "--runHiCpca", action="store_true",
    default=False)
parser.add_argument('-c', "--correlation", action="store_true",
    default=False)
args = vars(parser.parse_args())

## This script is created by Kevin Gau on 2018-05-02.
## In this script, I use HOMER to analyze Hi-C genome-wide
## interaction data from HiC-Pro output.

## The documentation of HOMER is at: 
## http://homer.ucsd.edu/homer/interactions/index.html

## our data:
##    cell type | cutter | replicate |   from    |  genome  |
## 1.   mESC    |  NCol  |     1     |  Ren Lab  |   mm9    |
## 2.   mESC    |  Mobi  |     1     |  Xie Lab  |   mm9    |
## 3.   mESC    |  Mobi  |     2     |  Xie Lab  |   mm9    |
## 4.   mESC    |  Mobi  |    WT1    |  Song Lab |   mm9    |
## 5.   mESC    |  Mobi  |    KO2    |  Song Lab |   mm9    |

HOMEPATH = "/home/data2/caojun/further_analysis"
#EXP_LIST = ['mESC_ren_NCol', 'mESC_song_DKO1', 'mESC_song_WT1', 
#'mESC_xie_500cellsRep1', 'mESC_xie_500cellsRep2']
EXP_LIST = ['PUM1-2','PUM2-2','DKO-2','WT-2']
RES_LIST = []

## makeTagDirectory - special paired-end operations for 
## making HOMER-style tag directories and filtering options
## for Hi-C
## Doc: 
## http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html
def makeTagDirectory(out_dir, in_file, res_site, genome):
    cmd = "makeTagDirectory %s -format HiCsummary %s \
-restrictionSite %s -genome %s" %(
        out_dir, in_file, res_site, genome)
    print cmd
    sys.stdout.flush()
    os.system(cmd)

if args["makeTagDirectory"]:
    print 'makeTagDirectory...'
    sys.stdout.flush()
    for exp in EXP_LIST:
        in_file = os.path.join(HOMEPATH, 'HiC-Pro_results/allValidPairs',
            exp+'_allValidPairs')
        out_dir = os.path.join(HOMEPATH, 'homerTag', exp)
        if 'mESC_ren_NCol' in in_file:
            makeTagDirectory(out_dir, in_file, 'CCATGG', 'mm9')
        else:
            makeTagDirectory(out_dir, in_file, 'GATC', 'mm9')


## analyzeHiC basic usage:
## analyzeHiC <Hi-C Tag Directory> [options] > outputMatrixFile.txt
## Correlation of normalized interaction ratios
## doc:
## http://homer.ucsd.edu/homer/interactions/HiCmatrices.html
resol = {
    '10k': 10000,
    '20k': 20000,
    '40k': 40000,
    '100k': 100000,
    '500k': 500000,
    '1M': 1000000
}

def getchrs(genome):
    if genome in ['mm9', 'mm10']:
        chrs = xrange(1,20)
        chrs = map(str, chrs)
        chrs.append('X')
    elif genome in ['hg19', 'hg38']:
        chrs = xrange(1,23)
        chrs = map(str, chrs)
        chrs.append('X')
    else:
        chrs = []
    return chrs

def calCorrMatrix(tag_dir, chrom, res, out_file):
    cmd = "analyzeHiC %s -chr chr%s -res %s -corr -cpu 8 >%s" %(
        tag_dir, chrom, res, out_file)
    print cmd
    sys.stdout.flush()
    os.system(cmd)

## runHiCpca basic usage:
## runHiCpca.pl <outputPrefix> <HiC Tag Directory> [options]
## doc:
## http://homer.ucsd.edu/homer/interactions/HiCpca.html
def runHiCpca(out_prefix, tag_dir, res, genome, active):
    cmd = "runHiCpca.pl %s %s -pc 2 -res %s -superRes %s -cpu 8 " %(
        out_prefix, tag_dir, res, res*2)
    if not active:
        active_op = "-genome %s" %genome
    else:
        active_op = "-active %s" %active
    cmd += active_op
    print cmd
    sys.stdout.flush()
    os.system(cmd)

if args["correlation"]:
    print "calculate correlation matrix..."
    sys.stdout.flush()
    for exp in EXP_LIST:
        tag_dir = os.path.join(HOMEPATH, 'homerTag', exp)
        for chrom in getchrs('mm9'):
            out_file = os.path.join(HOMEPATH, 'compartment/corrMatrix',
                exp+'_chr%s.corrMat100k.txt' %chrom)
            calCorrMatrix(tag_dir, chrom, resol['100k'], out_file)

if args['runHiCpca']:
    print "compartment analysis use PCA..."
    sys.stdout.flush()
    active = os.path.join(HOMEPATH, 'reference',
        'RenLab-H3K4me3-mESC-peaks.bed')
    for exp in EXP_LIST:
        tag_dir = os.path.join(HOMEPATH, 'homerTag', exp)
        out_dir = os.path.join(HOMEPATH, 'compartment/pca_100k')
        os.path.exists(out_dir) or os.makedirs(out_dir)
        out_prefix = os.path.join(out_dir, exp)
        runHiCpca(out_prefix, tag_dir, resol['100k'], 'mm9', active)


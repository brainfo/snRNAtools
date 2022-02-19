#!/home/hong/anaconda3/bin/python

## This script is for downloading fq files from json files ##
## I usually collected these json files from EBI website ##
## Creator: Hong Jiang ##
## 2021-10-21 at KI ##

import argparse
import json
from utils import parallel
import itertools
import pandas as pd

## Parameters
parser = argparse.ArgumentParser()

parser.add_argument('-n','--threads', action='store', dest='threads', type=int, help='threads')
parser.add_argument('-i','--input', action='store', dest='input', help='json or sdrf file')
parser.add_argument('-p','--folder', action='store', dest='folder', help='folder to store downloaded files')

paras = parser.parse_args()

def read_in_json(json_input):
    with open(json_input,'r') as json_to_load:
        info = json.load(json_to_load)
    fq_to_download = [srx['fastq_ftp'].split(';') for srx in info]
    return fq_to_download

def read_in_sdrf(sdrf_input):
    info = pd.read_csv(sdrf_input,sep='\t')
    fq_to_download = info.filter(regex='FASTQ_URI',axis=1).transpose().values.tolist()
    return fq_to_download

def download_list(fq_to_download):
    fqs = list(itertools.chain.from_iterable(fq_to_download))
    cmds = ['wget -nc -P {} '.format(paras.folder) + fq for fq in fqs]
    return cmds

if __name__ == '__main__':
    if 'json' in paras.input:
        cmds = download_list(read_in_json(paras.input))
    if 'sdrf' in paras.input:
        cmds = download_list(read_in_sdrf(paras.input))
    parallel.exe_parallel(cmds,paras.threads)

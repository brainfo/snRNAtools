import glob
import os
import pandas as pd
from utils import parallel
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()

parser.add_argument('-n','--threads', action='store', dest='threads', type=int, help='threads')
parser.add_argument('-i','--input_folder', action='store', dest='input_folder', help='input bam folder')
parser.add_argument('-o','--output_folder', action='store', dest='output_folder', help='output fragment folder')

paras = parser.parse_args()

def stat():
    cmds = []
    for bam in sorted(glob.glob(paras.input_folder+'/**nobl_clean_sort.bam')):
        name = '_'.join(bam.split('/')[-1].split('_')[:2])
        frag = f'{paras.output_folder}/{name}_frags.bed'
        cmd = [f'sinto fragments -b {bam} -f {frag} -m 10 -p 8 --use_chrom "^[0-9]+|X|Y" -t RG']
        cmds += cmd
    return cmds

if __name__ =='__main__':
    cmds = stat()
    parallel.exe_parallel(cmds, paras.threads)

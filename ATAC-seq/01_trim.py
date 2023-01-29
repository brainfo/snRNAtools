''' using NGmerge to trim adapters from paired-end reads '''

import glob
from utils import parallel
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-n','--threads', action='store', dest='threads', type=int, help='threads')
parser.add_argument('-i','--input_folder', action='store', dest='input_folder', help='input fastq folder')
parser.add_argument('-o','--output_folder', action='store', dest='output_folder', help='output fastq folder')

paras = parser.parse_args()

def trim():
    R1 = sorted(glob.glob(paras.input_folder+'/**/**_1.fq.gz'))
    cmds = []
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = r1.replace('_1.fq.gz', '_2.fq.gz', 1)
        name = '_'.join(r1.split('/')[-1].split('_')[:2])
        out1 = f'{paras.output_folder}/{name}_trim_1.fq.gz'
        out2 = f'{paras.output_folder}/{name}_trim_2.fq.gz'
        cmd = [f'bbduk.sh in1={r1} in2={r2} out1={out1} out2={out2} qtrim=rl trimq=20 maq=10 minlen=20 tbo ref=adapters.fa']
        cmds += cmd
    return cmds

if __name__ == '__main__':
    ## make output folder
    # os.makedirs(paras.output_folder)
    cmds = trim()
    parallel.exe_parallel(cmds,paras.threads)

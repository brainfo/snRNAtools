import glob
import os
import pandas as pd
from utils import parallel
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()

parser.add_argument('-n','--threads', action='store', dest='threads', type=int, help='threads')
parser.add_argument('-i','--input_folder', action='store', dest='input_folder', help='input fastq folder')
parser.add_argument('-o','--output_folder', action='store', dest='output_folder', help='output bam folder')

paras = parser.parse_args()

def get_CB():
    ## barcode dict
    bc_df = pd.read_csv('/mnt/data/hong/2022/mouse_oocyte/smartseq3/rawdata/reads_for_zUMIs.samples.txt', sep ='\t')
    bc_dict = defaultdict()
    for index, row in bc_df.iterrows():
        bc_dict[row['sample'].split('_')[1]] = row['BC']
    return bc_dict
    ## return barcode

def align():
    R1 = sorted(glob.glob(paras.input_folder+'/**_1.fq.gz'))
    cmds = []
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = r1.replace('_1.fq.gz', '_2.fq.gz', 1)
        name = '_'.join(r1.split('/')[-1].split('_')[:2])
        s_name = name.split('_')[1]
        cb_tag = bc_dict[s_name]
        out = f'{paras.output_folder}/{name}_sorted.bam'
        cmd = [f'bowtie2 --very-sensitive -k 10 -x /mnt/data/hong/reference/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/genome -1 {r1} -2 {r2} -p 8 --rg-id {name} --rg "CB:{cb_tag}"\
        | samtools view -u -@ 8\
        | samtools sort -o {out} -@ 8']
        cmds += cmd
    return cmds

def filter_bl():
    bl_file = '/mnt/data/hong/reference/Mus_musculus/NCBI/GRCm38/Sequence/mm10-blacklist.v2.bed'
    samtools_params = f'-M -L {bl_file} -O BAM -@ 8'
    cmds = []
    for bam in sorted(glob.glob(paras.output_folder+'/**_sorted.bam')):
        name = '_'.join(bam.split('/')[-1].split('_')[:2])
        nobl = f'{paras.output_folder}/{name}_nobl.bam'
        bl = f'{paras.output_folder}/{name}_bl.bam'
        # cmd = [f'samtools index {bam}; samtools view {samtools_params} --unoutput {nobl} {bam}']
        # cmd = [f'samtools view {samtools_params} --unoutput {nobl} -o {bl} {bam}']
        cmd = [f'bedtools intersect -abam {bam} -b {bl_file} -v > {nobl}']
        cmds += cmd
    return cmds

def filter_others():
    """
    use sambamba to filter out
    Reads that were not mapped, not primary alignment, missing a
    mate, mapq <10, or overlapping ENCODEs blacklist regions
    """
    # sambaster_params = 'mapping_quality >= 10 and not (unmapped or mate_is_unmapped) and first_of_pair and proper_pair'
    cmds = []
    for bam in sorted(glob.glob(paras.output_folder+'/**_nobl.bam')):
        name = '_'.join(bam.split('/')[-1].split('_')[:2])
        clean = f'{paras.output_folder}/{name}_nobl_clean.bam'
        clean_sort = f'{paras.output_folder}/{name}_nobl_clean_sort.bam'
        cmd = [f'sambamba view -h -t 8 -f bam \
            -F "paired and mapping_quality >= 10 and not (unmapped or mate_is_unmapped or secondary_alignment)" {bam} > {clean}; samtools sort -@ 8 {clean} > {clean_sort}; samtools index {clean_sort}']
        cmds += cmd
    return cmds

if __name__ =='__main__':
    #os.makedirs(paras.output_folder)
    # bc_dict = get_CB()
    # cmds_align = align()
    # parallel.exe_parallel(cmds_align, paras.threads)
    cmds_filter_bl = filter_bl()
    parallel.exe_parallel(cmds_filter_bl, paras.threads)
    cmds_filter_others = filter_others()
    parallel.exe_parallel(cmds_filter_others, paras.threads)
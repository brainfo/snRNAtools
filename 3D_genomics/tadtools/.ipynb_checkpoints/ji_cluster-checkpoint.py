import pandas as pd
import numpy as np
import operator
import sys
sys.path.insert(0, '/home/data4/jh/data/tools/customized/tadtools/tadtool_pkg')
import utils as ut
import itertools
import functools
import pickle

cluster_path = '/home/data4/jh/data/projects/2020/aggregate/mESC/E14TG2a/cluster/ipy/data/cluster6'

if __name__ == "__main__":
    print('reading file')
    npc = pd.read_csv('NPC_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    a11 = pd.read_csv(f'{cluster_path}/A00',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    a12 = pd.read_csv(f'{cluster_path}/A01_both',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    a22 = pd.read_csv(f'{cluster_path}/A11',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    s11 = pd.read_csv(f'{cluster_path}/B00',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    s12 = pd.read_csv(f'{cluster_path}/B01_both',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    s22 = pd.read_csv(f'{cluster_path}/B11',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
    oocyte = pd.read_csv('Oocyte_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    pn3 = pd.read_csv('PN3_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    sperm = pd.read_csv('Sperm_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    mesc_list = [a11,a12,a22,s11,s12,s22]
    for i, mesc in enumerate(mesc_list):
        tad_list = [mesc,npc,oocyte,sperm,pn3]
        tad_pair = list(itertools.combinations(tad_list,2))[:4]
        jis_list = list(map(lambda x: ut.compute_ji_ctr(x[0],x[1]),tad_pair))
        jis_file = open(f'./cluster6/mesc_{i}.obj','wb')
        pickle.dump(jis_list,jis_file)


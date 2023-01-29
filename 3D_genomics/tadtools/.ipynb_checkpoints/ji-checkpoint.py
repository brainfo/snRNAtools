import pandas as pd
import numpy as np
import operator
import sys
sys.path.insert(0, '/home/data4/jh/data/tools/customized/tadtools/tadtool_pkg')
import utils as ut
import itertools
import functools
import pickle

if __name__ == "__main__":
    print('reading file')
    npc = pd.read_csv('NPC_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    mesc = pd.read_csv('TADs.bed',sep='\t',header=None,names=['chrom','start','end'])
    oocyte = pd.read_csv('Oocyte_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    pn3 = pd.read_csv('PN3_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    sperm = pd.read_csv('Sperm_tads.bed',sep='\t',header=None,names=['chrom','start','end'])
    tad_list = [mesc,npc,oocyte,sperm,pn3]
    tad_pair = list(itertools.combinations(tad_list,2))
    print('start computing')
    jis_list = list(map(lambda x: ut.compute_ji_ctr(x[0],x[1]),tad_pair))
    jis_file = open('jis_test.obj','wb')
    print('start dumping')
    pickle.dump(jis_list,jis_file)
    

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
    print('start computing')
    jis_list = ut.compute_ji_ctr(mesc,npc)
    jis_file = open('test.obj','wb')
    print('start dumping')
    pickle.dump(jis_list,jis_file)
    

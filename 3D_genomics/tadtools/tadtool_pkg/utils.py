import numpy as np
import operator
import functools
from scipy.stats import mannwhitneyu
import itertools

chromosome = ['chr'+ i for i in [str(i) for i in range(1,20)]]
chromosome.append('chrX')

def out_boundary(insu,name,threshold):
    rep2Bound = insu.loc[insu.boundary_strength_1000000.notna(),]
    rep2BoundStrength = rep2Bound.loc[insu.boundary_strength_1000000>threshold,]
    rep2BoundStrength.to_csv(f'{name}_{threshold}_boundary.txt',sep='\t',index=False)
    return rep2BoundStrength

def filter_boundary(df):
    return df.loc[(df.chrom != 'chrY') & (df.is_bad_bin==False),:]

def to_tads(cluster,name,chromosome):
    data = cluster[['chrom','start','end']]
    for chrom in chromosome:
        data_chr = data.loc[data.iloc[:,0]==chrom,[False,True,True]].to_numpy().flatten()[1:-1]
        tad_chr = data_chr.reshape((-1,2))
        with open(f'{name}_tads.bed','a') as fo:
            for i in range(tad_chr.shape[0]):
                fo.write(f'{chrom}\t{tad_chr[i,0]}\t{tad_chr[i,1]}\n')

def set_intersect_len(set1,array):
    a, b = array
    intersect = set1.intersection(range(a, b))
    return len(intersect)
    
    
def compute_ji(tad1,tad2):
    jis = []
    for index, row in tad1.iterrows():
        tad1_tad = set(range(row[1],row[2]))
        tad2_chr = tad2.loc[tad2.chrom==row[0],['start','end']].to_numpy()
        intersect_lens = np.apply_along_axis(functools.partial(set_intersect_len,tad1_tad), 1, tad2_chr)
        if all(v == 0 for v in intersect_lens):
            jis.append(0)
        else:
            index, value = max(enumerate(intersect_lens), key=operator.itemgetter(1))
            union_len = len(tad1_tad.union(tad2_chr[index,:]))
            ji = value/union_len
            jis.append(ji)
    return jis 

def compute_ji_ctr(tad1,tad2):
    def body(element):
        index, row = element
        tad1_tad = set(range(row[1],row[2]))
        tad2_chr = tad2.loc[tad2.chrom==row[0],['start','end']].to_numpy()
        intersect_lens = np.apply_along_axis(functools.partial(set_intersect_len,tad1_tad), 1, tad2_chr)
        if all(v == 0 for v in intersect_lens):
            ji = 0
        else:
            index, value = max(enumerate(intersect_lens), key=operator.itemgetter(1))
            union_len = len(tad1_tad.union(tad2_chr[index,:]))
            ji = value/union_len
        return ji
    
    return list(map(body,tad1.iterrows()))

def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def sig(data,col_name):
    data_list = list(map(lambda x: x[col_name],data))
    data_pair = list(itertools.combinations(data_list,2))
    ps = []
    for pair in data_pair:
        stat, p = mannwhitneyu(pair[0], pair[1])
        ps.append(p)
    return ps

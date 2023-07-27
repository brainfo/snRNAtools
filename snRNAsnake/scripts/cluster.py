import pandas as pd
import scanpy as sc
import anndata as ad
import importlib
from sc_utils import scanpy_utils as su
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
from matplotlib.pyplot import rc_context
import numpy as np
import scipy.sparse as sp
import anndata
from typing import Dict, Optional
from collections import defaultdict
import itertools
import pygame
from matplotlib import cm
import matplotlib
import seaborn as sns
from collections import Counter

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_theme(style="ticks", rc=custom_params)
matplotlib.rcParams['figure.figsize'] = [6, 5]

sc.settings.verbosity = 0
sc.logging.print_header()
sc.set_figure_params(dpi=900)

## from snakemake input, output and config
workdir = snakemake.config["project"]["workdir"]
project_name = snakemake.config["project"]["project_name"]
h5dir = snakemake.config["project"]["h5dir"]

os.chdir(workdir)

filter_params = {
    'min_counts':400, 'min_genes':200, 'max_genes' : 5000, 'percent_mt':5, 'percent':3, 'filter_mt':True
}

raw_dict = defaultdict(lambda: "Not Present")
no_doublet_dict = defaultdict(lambda: "Not Present")
filter_dict = defaultdict(lambda: "Not Present")
## read in sample info
info = pd.read_csv(snakemake.input['info'], sep='\t')
cols = ['randnr_mother', 'randomization_mother']
group = info[cols].groupby('randomization_mother') \
                 .apply(lambda x: x.set_index('randomization_mother').to_dict('list')) \
                 .to_dict()

group_inv = defaultdict()
for k, v in group.items():
    for vi in v:
        group_inv[vi] = k
        
for root, sample_list, filenames in os.walk(f'{h5dir}'):
    for sample_name in sample_list:
        if sample_name.startswith('Placenta_'):
            raw_dict[sample_name] = su.anndata_from_h5(f'{h5dir}/{sample_name}/{sample_name}_cellbender_filtered.h5')

# doublet_params = {
#     'placenta_314':0.09, 'placenta_40':0.09, 'placenta_248' : 0.14, 'placenta_25':0.12, 
#     'placenta_226':0.12, 'placenta_60':0.15, 'placenta_373':0.11, 'placenta_32':0.12, 'placenta_303':0.12,
#     'placenta_357':0.14, 'placenta_312':0.13, 'placenta_306':0.12, 'placenta_330':0.13, 'placenta_81': 0.14
# }
for sample_name, sample in raw_dict.items():
    sample.var_names_make_unique('+')
    sc.external.pp.scrublet(sample)
    su.doublet_plot(sample_name, sample)

for sample_name, sample in raw_dict.items():
    doublet = np.array(sample.obs['predicted_doublet'], dtype=bool)
    no_doublet_dict[sample_name] = sample[~doublet]
for sample_name, sample in no_doublet_dict.items():
    su.qc(sample, f'{sample_name}_no_doublet', 'MT')
    filter_dict[sample_name] = su.filter_adata(sample, **filter_params)
ad_all = ad.concat(list(filter_dict.values()), label='sample', keys=list(filter_dict.keys()), join='outer', index_unique='-', merge='same')

## save ad_all
ad_all.write(snakemake.output['ad_all'], compression='gzip')
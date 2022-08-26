#!/mnt/data/hong/anaconda3/bin/python
## this is to merge fqs from one replicate ##
## but sequenced in different lanes ##
## created by Ruby Jiang on 2022-05-09 at KI ##
import os
from os import walk
fqs_folder = '../output/trim_glare/'
class replicate:
    r1s = []
    r2s = []
    names = {}
    def __init__(self, name, r1_suffix, r2_suffix):
        self.name = name
        self.r1_suffix = r1_suffix
        self.r2_suffix = r2_suffix
        self.r1_names = []
        self.r2_names = []
    def get_r1s(self):
        #print(f'{fqs_folder}{self.name}')
        for (dirpath, dirnames, filenames) in walk(f'{fqs_folder}{self.name}'):
            for filename in filenames:
                #print(filename, self.r1_suffix)
                if filename.endswith(self.r1_suffix):
                    #print(filename)
                    #self.r1s.append(filename)
                    self.r1_names.append(filename.removesuffix(self.r1_suffix))
        self.r1_names.sort()
        #print(self.r1_names)
    def get_r2s(self):
        for (dirpath, dirnames, filenames) in walk(f'{fqs_folder}{self.name}'):
            for filename in filenames:
                if filename.endswith(self.r2_suffix):
                    #self.r2s.append(filename)
                    self.r2_names.append(filename.removesuffix(self.r2_suffix))
        self.r2_names.sort()
    def get_names(self):
        print(self.r1_names, self.r2_names)
        assert self.r1_names == self.r2_names, "paired reads do not match"
        self.names = set(self.r1_names)
        assert len(self.names) > 0, "No sequence file found"
    def concate_fqs(self):
        self.r1s = [fqs_folder+self.name+'/'+r+self.r1_suffix for r in self.names]
        self.r2s = [fqs_folder+self.name+'/'+r+self.r2_suffix for r in self.names]
        if len(self.r1s) > 1:
            r1_zcat_list = ' '.join(self.r1s)
            cmd1 = f'cat {r1_zcat_list} > {self.name}_1.fq.gz'
            os.system(cmd1)
            r2_zcat_list = ' '.join(self.r2s)
            cmd2 = f'cat {r2_zcat_list} > {self.name}_2.fq.gz'
            os.system(cmd2)
## could write as get jobs and then submit in parellel
for (dirpath, dirnames, filenames) in walk(fqs_folder):
    for dirname in dirnames:
        pgc_rep = replicate(dirname, '_1_val_1.fq.gz', '_2_val_2.fq.gz')
        pgc_rep.get_r1s()
        pgc_rep.get_r2s()
        pgc_rep.get_names()
        pgc_rep.concate_fqs()





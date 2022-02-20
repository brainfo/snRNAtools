''' using snpsplit to seperate allels '''

import glob
import os
import sys

cwd = os.getcwd()

def outFolder():
    inputpath = os.path.join(cwd,'allelbam')
    outputpath = os.path.join(cwd,'snpsplit')
    print(outputpath)

    for dirpath, dirnames, filenames in os.walk(inputpath):
        print(dirpath)
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            print(structure)
            os.mkdir(structure)
        else:
            print("Folder does already exits!")

def snp():
    bamin = sorted(glob.glob(cwd+'/allelbam/*/*/*_noMT.bam'))
    
    for i in range(len(bamin)):
        bi = bamin[i]
        bo = bi.replace('allelbam','snpsplit')
        cmd = 'samtools sort -n -o %s-@ 12 %s' %(bo,bi)
        os.system(cmd)
        dupbam = '/'.join(bo.split('/')[:-1]) + '/nodup.'+bo.split('/')[-1]
        cmd = 'sambamba view -h -t 12 -f bam -F "not ([XA] != null or [SA] != null) and not unmapped and not duplicate" %s -o %s' %(bo,dupbam)
        os.system(cmd)
        filbam = '/'.join(bo.split('/')[:-1]) + '/filter.'+bo.split('/')[-1]
        cmd = 'bedtools intersect -v -abam %s -b %s > %s' %(dupbam,black,filbam)
        os.system(cmd)
        cmd = 'samtools index %s' %(filbam)
        os.system(cmd)
        MTbam = '.'.join(filbam.split('.')[:-1])+'_noMT.bam'
        cmd = 'samtools idxstats %s | cut -f 1 | grep -v MT | xargs samtools view -b %s > %s' %(filbam,filbam,MTbam)
        os.system(cmd)
        cmd = 'samtools index %s' %(MTbam)
        os.system(cmd)

if __name__ == '__main__':
    outFolder()
    align()
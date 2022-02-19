''' using NGmerge to trim adapters from paired-end reads '''

import glob
import os
import sys

cwd = os.getcwd()

def outFolder():
    inputpath = os.path.join(cwd,'fq')
    outputpath = os.path.join(cwd,'bbduk')

    for dirpath, dirnames, filenames in os.walk(inputpath):
        structure = os.path.join(outputpath, dirpath[len(inputpath):])
        if not os.path.isdir(structure):
            os.mkdir(structure)
        else:
            print("Folder does already exits!")


def trim():
    R1 = sorted(glob.glob(cwd+'/fq/**/**/*_1.fastq.gz'))
    R2 = sorted(glob.glob(cwd+'/fq/**/**/*_1.fastq.gz'))
    
    assert len(R1)==len(R2)
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = R2[i]
        out1 = r1.replace('fq','bbduk',1)
        out2 = r2.replace('fq','bbduk',1)
        log = path + 'log/'+nameR1+'.log'
        cmd = 'bbduk.sh in1=%s in2=%s out1=%s out2=%s qtrim= trimq=10 maq=10 minlen=20' %(r1,r2,out1,out2)
        os.system(cmd)

if __name__ == '__main__':
    outFolder()
    trim()

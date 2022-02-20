''' using NGmerge to trim adapters from paired-end reads '''

import glob
import os
import sys

cwd = os.getcwd()

def outFolder():
    inputpath = os.path.join(cwd,'fq')
    outputpath = os.path.join(cwd,'bbduk')
    print(outputpath)

    for dirpath, dirnames, filenames in os.walk(inputpath):
        print(dirpath)
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            print(structure)
            os.mkdir(structure)
        else:
            print("Folder does already exits!")


def trim():
    R1 = sorted(glob.glob(cwd+'/fq/2cell/**/*_1.fastq.gz'))
    R2 = sorted(glob.glob(cwd+'/fq/2cell/**/*_2.fastq.gz'))
    assert len(R1)==len(R2)
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = R2[i]
        #print(r1)  
        out1 = r1.replace('fq','bbduk',1)
        out2 = r2.replace('fq','bbduk',1)
        name = '_'.join(r1.split('/')[-1].split('_')[:-1])
        #print(out1)
        log = os.path.join(cwd,'log') + '/' + name + '_bbduk.log'
        err = os.path.join(cwd,'log') + '/' + name + '_bbduk.err'
        cmd = 'nohup bbduk.sh in1=%s in2=%s out1=%s out2=%s qtrim=r trimq=10 maq=10 minlen=20 >%s 2>%s &' %(r1,r2,out1,out2,log,err)
        os.system(cmd)

if __name__ == '__main__':
    #outFolder()
    trim()

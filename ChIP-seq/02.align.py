''' using bwa to align paired ends '''

import glob
import os
import sys

cwd = os.getcwd()
ref = '/home/data4/jh/data/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/maskIndex/mm10_C57bl6_DBA2J'
ref_bwa = '/home/data4/jh/data/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'
def outFolder(fo):
    inputpath = os.path.join(cwd,'bbduk')
    outputpath = os.path.join(cwd,fo)
    print(outputpath)

    for dirpath, dirnames, filenames in os.walk(inputpath):
        print(dirpath)
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            print(structure)
            os.mkdir(structure)
        else:
            print("Folder does already exits!")


def align():
    R1 = sorted(glob.glob(cwd+'/bbduk/8cell/**/**/*_1.fastq.gz'))
    R2 = sorted(glob.glob(cwd+'/bbduk/8cell/**/**/*_2.fastq.gz'))
    
    assert len(R1)==len(R2)
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = R2[i]
        outPath = '/'.join(r1.replace('bbduk','bwa').split('/')[:-1])
        name = '_'.join(r1.split('/')[-1].split('_')[:-1]) + '.bam'
        #print(name)
        out = os.path.join(outPath,name)
        cmd = 'bwa mem %s %s %s -t 12 | samtools view -h -S -b -o %s -@ 12 2>./log/%s.err' %(ref_bwa,r1,r2,out,name)
        os.system(cmd)

def align_bowtie2():
    R1 = sorted(glob.glob(cwd+'/bbduk/2cell/**/*_1.fastq.gz'))
    R2 = sorted(glob.glob(cwd+'/bbduk/2cell/**/*_2.fastq.gz'))
    
    assert len(R1)==len(R2)
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = R2[i]
        outPath = '/'.join(r1.replace('bbduk','bowtie2').split('/')[:-1])
        name = 'bowtie_'+'_'.join(r1.split('/')[-1].split('_')[:-1]) + '.bam'
        print(name)
        out = os.path.join(outPath,name)
        cmd = 'bowtie2 -p 14 -q --local \
        -x %s -1 %s -2 %s -S | \
        samtools view -h -S -b -o %s' %(ref,r1,r2,out)
        os.system(cmd)


if __name__ == '__main__':
    #outFolder('bwa')
    align()

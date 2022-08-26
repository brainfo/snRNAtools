import glob
import os
import sys

path = os.getcwd()
in_path = path +'/NGmerge/'
print(path)

def align(path):
    R1 = sorted(glob.glob(in_path+'/*_1.fastq.gz'))
    R2 = sorted(glob.glob(in_path+'/*_2.fastq.gz'))
    assert len(R1)==len(R2)
    for i in range(len(R1)):
        r1 = R1[i]
        r2 = R2[i]
        nameR1 = r1.split('/')[-1].split('_')[0]
        nameR2 = r2.split('/')[-1].split('_')[0]
        assert nameR1==nameR2
        print(nameR1)
        cmd = 'bowtie2 --very-sensitive -k 10 -x /home/data4/jh/data/reference/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -1 %s -2 %s -p 8\
        | samtools view -u -@ 8\
        | samtools sort -n -o %s -@ 8' %(r1,r2,path+'/align/'+nameR1+'.bam')
        os.system(cmd)

if __name__ =='__main__':
    align(path)
import glob
import os,sys,time


path = os.getcwd()
in_path = path +'/result/align'

def Dsample(sample):
    flag = './bam/%s.bam.flagstat' %sample
    bam = './bam/%s.bam' %sample
    with open(flag,'r') as f:
        line = f.readlines()[4]
        match_number = re.match(r'(\d.+) \+.+',line)
        total_reads = int(match_number.group(1))

    target_reads = 1455932
    if total_reads >= target_reads:
        down_rate = target_reads/total_reads
    else:
        down_rate = 1

    cmd = 'sambamba view -f bam -t 5 --subsampling-seed=3 -s %s %s | samtools sort -@ 5 -T %s.tmp > ./bam/subsample.%s.bam 2> ./log/%s.subsample.log' %(down_rate,bam,sample,sample,sample)
    os.system(cmd)
    cmd = 'nohup samtools index ./bam/subsample.%s.bam &' %sample
    os.system(cmd)

def makeBW(sample):
    nameR1 = r1.split('/')[-1].split('.')[0]
    cmd = 'bamCoverage -b %s --normalizeUsing RPGC--binSize 30 --smoothLength 300 -p 16 --extendReads 200 --effectiveGenomeSize 2913022398 --extendReads --ignoreDuplicates --ignoreForNormalization chrX -o ./bw/%s.bw 2> ./log/%s.bw.log' %()
    os.system(cmd)

if __name__ == '__main__':
    Rep1 = sorted(glob.glob(in_path+'/*.bam'))
    for i in range(len(Rep1)):
        r1 = Rep1[i]
        makeBW(r1)

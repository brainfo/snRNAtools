''' call peaks from bam files using macs2 '''

import glob
import os
import sys

cwd = os.getcwd()

def outFolder():
    inputpath = os.path.join(cwd,'bam')
    outputpath = os.path.join(cwd,'peak')
    print(outputpath)

    for dirpath, dirnames, filenames in os.walk(inputpath):
        print(dirpath)
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            print(structure)
            os.mkdir(structure)
        else:
            print("Folder does already exits!")


def macs2():
    bamin = sorted(glob.glob(cwd+'/bam/**/**/**/*.bam'))
    for bam in bamin:
        if 'input' not in bam:
            if 'merge' in bam and 'nodup' in bam:
                ctlpath = '/'.join(bam.split('/')[:-2])+'/input'
                ctl = glob.glob(os.path.join(ctlpath,'nodup.*.bam'))[0]
                name = bam.split('/')[-1].strip('.bam')
                peakpath = '/'.join(bam.split('/')[:-1]).replace('bam','peak',1)
                cmd = 'macs2 callpeak -t %s -c %s -f BAM -g 2652783500 -n %s --outdir %s 2> ./log/macs2.%s.log' %(bam,ctl,name,peakpath,name)
                os.system(cmd)
            if 'merge' not in bam and 'filter' in bam:
                ctlpath = '/'.join(bam.split('/')[:-2])+'/input'
                ctl = glob.glob(os.path.join(ctlpath,'filter.*.bam'))[0]
                name = bam.split('/')[-1].strip('.bam')
                peakpath = '/'.join(bam.split('/')[:-1]).replace('bam','peak',1)
                cmd = 'macs2 callpeak -t %s -c %s -f BAM -g 2652783500 -n %s --outdir %s 2> ./log/macs2.%s.log' %(bam,ctl,name,peakpath,name)
                os.system(cmd)               


if __name__ == '__main__':
    #outFolder()
    macs2()


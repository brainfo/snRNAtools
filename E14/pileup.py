## calculate ChIP-seq pileup track ##
## Edited by Ruby Jiang ##
## This is very fast ##

import re
import os
import sys
import glob

cwd = os.getcwd()

def outFolder():
    inputpath = os.path.join(cwd,'bams')
    outputpath = os.path.join(cwd,'coverage')
    print(outputpath)

    for dirpath, dirnames, filenames in os.walk(inputpath):
        print(dirpath)
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            print(structure)
            os.mkdir(structure)
        else:
            print("Folder does already exits!")

def predictD(bam,name):
	dfile = './log/' + name + '_d.log'
	cmd = f'macs2 predictd -i {bam} -g mm --outdir ./log --rfile {name} 2> {dfile}'
	os.system(cmd)
	return dfile

def trieve_d(fi):
	with open(fi,'r') as file:
		data = file.read().replace('\n','')
		result = re.findall(r'predicted fragment length is (\d+) bps',data)
	return result
	##return a list of string

def pileup(bam,name):
	infile = bam
	#name = bam.split('/')[-1].strip('.bam')
	outfile = f'./coverage/{name}_RPKM.bedgraph'
	## macs2 dedup result is bdg format
	dfile = predictD(bam,name)
	d = trieve_d(dfile)[0]
	log = './log/' + name + '_coverage.log'
	#cmd = 'nohup macs2 pileup -i %s -o %s --extsize %s -f BEDPE 2>%s &' %(infile,outfile,d,log)
	## Not for paired-end
	cmd = 'bamCoverage -b %s -o %s -of bedgraph -p 16 --normalizeUsing RPKM -e %s --ignoreDuplicates 2>%s' %(infile,outfile,d,log)
	os.system(cmd)

if __name__ == '__main__':
	#outFolder()
	bams = sorted(glob.glob(cwd+'/bams/*.bam'))
	## this format is bedpe actually.
	for bam in bams:
		name = bam.split('/')[-1].strip('.bam')
		pileup(bam,name)

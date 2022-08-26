#!/home/hong/anaconda3/bin/python

## This util is for support omics analysis ##
## Especially for WGBS ##
## Creator: Hong Jiang ##
## 2022-03-08 at KI ##

import os

def outFolder(inputpath, outputpath):
    '''
    create same output folder structure as input
    '''
    for dirpath, dirnames, filenames in os.walk(inputpath):
        structure = os.path.join(outputpath, dirpath[len(inputpath)+1:])
        if not os.path.isdir(structure):
            os.mkdir(structure)
        else:
            print("Folder does already exits!")
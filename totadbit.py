
with open('RenLab-HiC-mESC-HindIII.allValidPairs', 'r') as fi:
    with open('test.txt', 'a') as fo:
        for line in fi:
            l = line.strip().split('\t')
            o = l[0:4]
            length = str(int(l[7])/2)
            o = o + [length,'0','0'] + l[7:10] + [length,'0','0']
            fo.write('\t'.join(o))
            fo.write('\n')
            

import glob
import os,sys,time
import threading

path = os.getcwd()
in_path = path +'/result/align'
print(path)

threads = 16

class execmd(threading.Thread):
    def __init__(self, cmd):
        super(execmd, self).__init__()
        self.cmd = cmd
    def run(self):
        now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        print('%s %s' %(now, self.cmd))
        os.system(self.cmd)

# def get_task_batch(path):
#     task_list = []
#     i,j = 0,0
#     Rep1 = sorted(glob.glob(in_path+'/*_rep1.bam'))
#     Rep2 = sorted(glob.glob(in_path+'/*_rep2.bam'))
#     #Rep3 = sorted(glob.glob(in_path+'/*_rep3.bam'))
#     for i in range(len(Rep1)):
#         #r3 = []
#         r1 = Rep1[i]
#         r2 = Rep2[i]
#         nameR1 = '_'.join(r1.split('/')[-1].split('_')[:-1])
#         nameR2 = '_'.join(r2.split('/')[-1].split('_')[:-1])
#         assert nameR1==nameR2
#         #r3 = [s for s in Rep3 if nameR1 in s]
#         #if len(r3)>0:
#             #print(nameR1,' 3 replicates')
#             #cmd = 'Genrich -t %s,%s,%s -o %s -j -y -r -e chrM -v -f %s' %(r1,r2,r3[0],path+'/result/peak/'+nameR1+'.narrowPeak',path+'/result/peak/'+nameR1+'.log')
#         #else:
#         print(nameR1,' 2 replicates')
#         cmd = 'Genrich -t %s,%s -o %s -j -y -r -e chrM -v -f %s' %(r1,r2,path+'/result/peak1/'+nameR1+'.narrowPeak',path+'/result/peak/'+nameR1+'.log')
#         print(j)

#         if j ==0:
#             task_list.append([])
#         i = j // threads
#         task_list[i].append(cmd)
#         j += 1
#     return task_list

##  no rep
def get_task_batch(path):
    task_list = []
    i,j = 0,0
    Rep1 = sorted(glob.glob(in_path+'/*.bam'))
    #Rep3 = sorted(glob.glob(in_path+'/*_rep3.bam'))
    for i in range(len(Rep1)):
        #r3 = []
        r1 = Rep1[i]

        nameR1 = r1.split('/')[-1].split('.')[0]
        #r3 = [s for s in Rep3 if nameR1 in s]
        #if len(r3)>0:
            #print(nameR1,' 3 replicates')
            #cmd = 'Genrich -t %s,%s,%s -o %s -j -y -r -e chrM -v -f %s' %(r1,r2,r3[0],path+'/result/peak/'+nameR1+'.narrowPeak',path+'/result/peak/'+nameR1+'.log')
        #else:
        print(nameR1,' 2 replicates')
        cmd = 'Genrich -t %s -o %s -j -y -r -e chrM -v -f %s' %(r1,path+'/result/peak/'+nameR1+'.narrowPeak',path+'/result/peak/'+nameR1+'.log')
        print(j)

        if j ==0:
            task_list.append([])
        i = j // threads
        task_list[i].append(cmd)
        j += 1
    return task_list
 
if __name__ == '__main__':
    for task_list in get_task_batch(path):
        for cmd in task_list:
            t = execmd(cmd)
            t.start()
        t.join()
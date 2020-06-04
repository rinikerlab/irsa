import numpy as np
import os
import pickle
import argparse
"""
Convert ORCA files to input files for the IRSA algorithm
"""
parser = argparse.ArgumentParser(description='Input converter')
parser.add_argument('-d', '--dir',required=False,type=str,default="./")
parser.add_argument('-n', '--nr',required=True,type=int)
args = parser.parse_args()
args = vars(parser.parse_args())

nr_at = args['nr']
directory = args['dir']

freq_big = []
energy_big = []
gibbs_big = []
files = []
for filename in os.listdir(directory):
    if filename.endswith(".engrad"): 
        f_open = open(directory+filename,"r")
        for i in range(7):
            f_open.readline()
        E = float(f_open.readline())
        f_open.close()
    if filename.endswith(".hess"): 
        f_open = open(directory+filename,"r")
        freq = [[],[]]
        ir=False
        for line in f_open:
            if("Final Gibbs free energy" in line):
                Gibbs = float(line.split()[-2])
            elif("IR SPECTRUM" in line):
                ir=True
                counter_tmp = 0
                f_open.readline()
                f_open.readline()
                f_open.readline()
                f_open.readline()
                while(counter_tmp<nr_at*3-6-5):
                    tmp = f_open.readline().split()
                    freq[0].append((float)(tmp[1]))
                    freq[1].append((float)(tmp[2]))
                    if(float(tmp[1])<0):
                        ir=False
                    counter_tmp+=1
        if(ir==False):
            cwd = os.getcwd()
        else:
            freq = np.swapaxes(np.asarray(freq),0,1)
            freq_big.append(np.asarray(freq))
            energy_big.append(np.asarray(E))
            gibbs_big.append(np.asarray(Gibbs))
            files.append(filename)
        f_open.close()

freq_big = np.asarray(freq_big).reshape(-1,nr_at*3-6-5,2)
energy_big = np.asarray(energy_big)
gibbs_big = np.asarray(gibbs_big)
print(energy_big.shape)

pickle.dump(energy_big,open("energy.p","wb"))
pickle.dump(freq_big,open("freq.p","wb"))
pickle.dump(gibbs_big,open("gibbs.p","wb"))
files = np.asarray(files)
idx = np.argsort(energy_big)

energy_big = energy_big[idx]*2625.50
energy_big -= np.min(energy_big)
gibbs_big = gibbs_big[idx]*2625.50
gibbs_big -= np.min(gibbs_big[idx])
files = files[idx]
f = open("sorted","w+")
for i in range(0,len(files)):
    f.write(str(energy_big[i])+" "+str(gibbs_big[i])+" "+files[i]+"\n")
f.close()

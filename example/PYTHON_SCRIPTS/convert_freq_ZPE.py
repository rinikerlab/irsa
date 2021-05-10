import numpy as np
import os
import pickle
freq_big = []
energy_big = []
gibbs_big = []
name = []
EE = 0
ZPE = 0
TVC = 0
TRC = 0
TTC = 0
FET = 0
kBT =  0.00094421
corr = 1.0192
for filename in os.listdir("./"):
    if filename.endswith(".out"): 
        f_open = open(filename,"r")
        freq = [[],[]]
        ir=False
        for line in f_open:
            if("FINAL SINGLE POINT ENERGY" in line):
                E = float(line.split()[-1])
            elif("Electronic energy" in line):
                EE = (float)(line.split()[3])
            elif("Zero point energy" in line):
                ZPE = (float)(line.split()[4])#*1.0187
            elif("Thermal vibrational correction" in line):
                TVC = (float)(line.split()[4])
            elif("Thermal rotational correction" in line):
                TRC = (float)(line.split()[4])
            elif("Thermal translational correction" in line):
                TTC = (float)(line.split()[4])
            elif("Total entropy correction" in line):
                FET = (float)(line.split()[4])
            elif("IR SPECTRUM" in line):
                ir=True
                counter_tmp = 0
                name.append(filename)
                f_open.readline()
                f_open.readline()
                f_open.readline()
                f_open.readline()
                while(True):
                    tmp = f_open.readline().split()
                    if(len(tmp)>=6):
                        print(tmp[1])
                        freq[0].append((float)(tmp[1]))
                        freq[1].append((float)(tmp[2]))
                        if(float(tmp[1])<0):
                            ir=False
                        counter_tmp+=1
                    else:
                        break
        if(ir==False):
            cwd = os.getcwd()
            print(filename)
        else:
            tmp = np.zeros((600,2))
            freq = np.swapaxes(np.asarray(freq),0,1)
            tmp[:len(freq),:] = freq[:,:]
            freq_big.append(np.asarray(tmp))
            energy_big.append(np.asarray(E))
            Gibbs = ZPE*corr+TRC+TTC+TVC+FET+kBT+E
            gibbs_big.append(np.asarray(Gibbs))
idx = np.argsort(np.asarray(gibbs_big))

freq_big = np.asarray(freq_big)#.reshape(-1,nr_at*3-6-5,2)
energy_big = np.asarray(energy_big)
gibbs_big = np.asarray(gibbs_big)
name = np.asarray(name)
pickle.dump(energy_big,open("energy.p","wb"))
pickle.dump(freq_big,open("freq.p","wb"))
pickle.dump(gibbs_big,open("gibbs.p","wb"))

import numpy as np
import os
import pickle
import argparse
"""
Convert GAUSSIAN09 files to input files for the IRSA algorithm
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
vcd_big = []
nr_at = 0
for filename in os.listdir(directory):
    if filename.endswith(".log"): 
        f_open = open(directory+"/"+filename,"r")
        freq = [[],[]]
        vcd = [[],[]]
        nr_at = 0
        v_tmp = []
        I_tmp = []
        R_tmp = []
        for line in f_open:
            if("SCF Done:" in line):
                E = float(line.split()[4])
            elif("NAtoms=" in line):
                nr_at=(int)(line.split()[1])
            elif("Frequencies --" in line):
                tmp = line.split()[2:]
                for el in tmp:
                    v_tmp.append(float(el))
                    if(float(el) < 0):
                        print("NEGATIVE FREQUENCIES")
                        print(os.getcwd())
            elif("IR Inten" in line):
                tmp = line.split()[3:]
                for el in tmp:
                    I_tmp.append(float(el))
            elif("Rot. str." in line):
                tmp = line.split()[3:]
                for el in tmp:
                    R_tmp.append(float(el))
        freq = np.zeros((len(I_tmp),2))
        vcd = np.zeros((len(I_tmp),2))
        freq[:,0] = np.asarray(v_tmp)
        freq[:,1] = np.asarray(I_tmp)
        vcd[:,0] = np.asarray(v_tmp)
        vcd[:,1] = np.asarray(R_tmp)
        freq_big.append(np.asarray(freq))
        energy_big.append(np.asarray(E))
        vcd_big.append(np.asarray(vcd))
        #gibbs_big.append(np.asarray(Gibbs))
freq_big = np.asarray(freq_big)#.reshape(-1,nr_at*3-6,2)
energy_big = np.asarray(energy_big)
vcd_big = np.asarray(vcd_big)
pickle.dump(energy_big,open("energy.p","wb"))
pickle.dump(freq_big,open("freq.p","wb"))
pickle.dump(vcd_big,open("vcd.p","wb"))

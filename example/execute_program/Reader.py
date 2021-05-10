import numpy as np
import pickle 

"""
This class reads in the data necessary for the algorithm.

The theoretical frequency spectra need to be dumped to a 
pickle file and must have the dimension (n_conformer,n_freq,2),
    where the first dimension is the number of conformer
    the second dimension is the number of normal modes, and
    the third dimension is the frequency (0) and dipol moment (1) OR ir intensity.
It must be specified, whether the IR intensity is saved or the dipole moment (default is IR intensity)

The energies of the conformer also must be dumped to a pickle file and must have the dimension (n_conformer).
The default energy assumed is Hartree (default in most QM programs), and is converted to kJ/mol.

"""
class Reader:
    def __init__(self):
        print("Reader init successfully")
    def read_pickle(self,file_name_E,file_name_f,threshold = 20):
        self.E = np.asarray(pickle.load(open(file_name_E,'rb')))
        self.freq = np.asarray(pickle.load(open(file_name_f,'rb')))
        ## SORT ENERGIES AND FREQUENCIES
        idx = np.argsort(self.E)
        self.E = self.E[idx]
        self.E -= np.min(self.E)
        self.freq = self.freq[idx]
        ## CONVERT HARTREE TO KJ/MOL
        self.E *= 2625.50
        print(self.E)
        tmp = self.E# threshold
        ##Only consider contributing conformers
        #self.E = self.E[0:1]
        #self.freq = self.freq[0:1]
    def get_freq(self):
        return self.freq
    def get_energy(self):
        return self.E 
    def write(self,out,p,sp,returnvalue,args):
        f_open = open(out+"results.txt","w+")
        f_open.write("PEARSON "+str(p)+"\n")
        f_open.write("SPEARMAN "+str(sp)+"\n")
        f_open.write("ALIGNMENT_SCORE "+str(-returnvalue)+"\n")
        f_open.write("PRODUCT "+str(-returnvalue*p)+"\n")
        f_open.write("------------SETTINGS USED---------------\n")
        f_open.write(str(args))
        f_open.close()

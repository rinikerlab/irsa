import sys
import random
import numpy as np
import argparse
import Reader
import Spectrum
import Algorithm

parser = argparse.ArgumentParser(description='Align Spectra')
parser.add_argument('-f', '--freq',required=True)
parser.add_argument('-e', '--energy',required=True)
parser.add_argument('-exp', '--experimental',required=True)

parser.add_argument('-w1','--width1',required=False,default=16,type=int)
parser.add_argument('-w2','--width2',required=False,default=16,type=int)
parser.add_argument('-d','--dipol',required=False,default=0,type=int)
parser.add_argument('-m','--mu',required=False,default=1.0,type=float)
parser.add_argument('-s0','--sigma0',required=False,default=0.3,type=float)
parser.add_argument('-s1','--sigma1',required=False,default=0.01,type=float)
parser.add_argument('-c','--cutoff',required=False,default=20,type=float)
parser.add_argument('-T','--temperature',required=False,default=298.15,type=float)
parser.add_argument('-lb','--lower_bound',required=False,default=1000,type=int)
parser.add_argument('-hb','--higher_bound',required=False,default=1800,type=int)
parser.add_argument('-res','--resolution',required=False,default=12,type=int)
parser.add_argument('-d0','--dummy0',required=False,default=1,type=int)
parser.add_argument('-d1','--dummy1',required=False,default=0,type=int)
parser.add_argument('-o','--out',required=False,default="",type=str)
args = parser.parse_args()
args = vars(parser.parse_args())
file_path_frequency = args['freq']
file_path_energy = args['energy']
file_path_exp = args['experimental']
width1 = args['width1']
width2 = args['width2']
dipol = args['dipol']
mu = args['mu']
sigma0 = args['sigma0']
sigma1 = args['sigma1']
cutoff = args['cutoff']
T = args['temperature']
resolution = args['resolution']
u = args['lower_bound']
h = args['higher_bound']
dummy_0 = args['dummy0']
dummy_1 = args['dummy1']
out = args['out']
#file_path_exp = args['-exp']
class Main:
    def __init__(self):
        print("ALIGNMENT STARTED")
    def run(self):
            self.Reader = Reader.Reader()
            self.Reader.read_pickle(file_path_energy,file_path_frequency)
            self.exp = np.loadtxt(file_path_exp,usecols=(0,1,))
            idx = np.argsort(self.exp[:,0])
            self.exp = np.asarray(self.exp[idx])
            self.freq = self.Reader.get_freq()
            self.E = self.Reader.get_energy()
            self.Spectrum = Spectrum.Spectrum(self.exp,self.freq,self.E,T=T,lower_bound=u,higher_bound=h,w1=width1,w2=width2,resolution_exp=resolution,Dipol=dipol,out=out)
            Alg = Algorithm.Algorithm(self.Spectrum,mu=mu,sigma_0=sigma0,sigma_1=sigma1,cutoff=cutoff,dummy_0=dummy_0,dummy_1=dummy_1)
            returnvalue, freq_old, freq_new, inten_new = Alg.Needleman_IR()
            print(-returnvalue)
            self.Spectrum.set_new_freq(freq_old,freq_new,inten_new)
            self.Spectrum.plot_unshifted()
            self.Spectrum.create_shifted()
            self.Spectrum.plot()
            p = self.Spectrum.integrate()
            self.Reader.write(out,p,returnvalue,args)
            return 1
        #except:
        #    return 1
program = Main()
if(program.run()):
    exit()
else:
    print("ALIGNMENT FAILED")


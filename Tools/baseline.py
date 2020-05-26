import numpy as np
import scipy
import matplotlib.pyplot as py
import argparse
"""
Perform a baseline correction for the experimental IR spectrum
"""
parser = argparse.ArgumentParser(description='Baseline correction')
parser.add_argument('-exp', '--experimental',required=True,type=str)
parser.add_argument('-l','--lower',required=False,type=int,default=1000)
parser.add_argument('-h','--higher',required=False,type=int,default=1600)

parser.add_argument('-d','--degree',required=False,type=int,default=0)
parser.add_argument('-o','--out',required=False,type=str,default="baseline.txt")

args = parser.parse_args()
args = vars(parser.parse_args())

exp = args['exp']
lb = args['lower']
hb = args['higher']
degree = args['degree']
out = args['out']

exp = np.loadtxt(exp,usecols=(0,1,))
vals = (exp[:,0] >= lb) & (hb >=exp[:,0])

p = np.polyfit(exp[vals,0],exp[vals,1]-np.min(exp[vals,1]),degree)
p = np.poly1d(p)
arr = p(exp[vals,0])
exp[vals,1] -= arr
exp[vals,1] -= np.min(a[vals,1])
py.plot(exp[vals,0],exp[vals,1],"baseline_corrected")
py.plot(exp[vals,0],arr,"baseline")
py.legend()
py.show()
np.savetxt(out,exp[vals])

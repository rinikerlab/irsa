#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
ir=np.loadtxt('IR_theo.txt',usecols=(0,1,))
idx = (ir[:,0]<2000) & (ir[:,0]>500)
ir=ir[idx]
ir[:,1]=ir[:,1]/np.max(ir[:,1])
idx = np.argsort(ir[:,0])
ir = ir[idx]
#idx = np.argsort(peaks[:,0])
xdat = ir[:,0]
ydat = ir[:,1]
peaks = []
peaks_y = []
for i in range(1,len(ir)-1):
    if(ir[i-1,1]<ir[i,1]>ir[i+1,1]): #and ir[i,0]<2000):
        peaks.append(ir[i,0])
        peaks_y.append(ir[i,1])
f=open('peaks_mod_theo.txt','w')
for i in range(len(peaks)):
    f.write(str(peaks_y[i])+" "+str(peaks[i])+" 0 \n")
f.close()


ir=np.loadtxt('IR.txt',usecols=(0,1,))
idx = (ir[:,0]<2000) & (ir[:,0]>500)
ir=ir[idx]
ir[:,1]=ir[:,1]/np.max(ir[:,1])
idx = np.argsort(ir[:,0])
ir = ir[idx]
#idx = np.argsort(peaks[:,0])
xdat = ir[:,0]
ydat = ir[:,1]
peaks = []
peaks_y = []
for i in range(1,len(ir)-1):
    if(ir[i-1,1]<ir[i,1]>ir[i+1,1]): #and ir[i,0]<2000):
        peaks.append(ir[i,0])
        peaks_y.append(ir[i,1])
f=open('peaks_mod_exp.txt','w')
for i in range(len(peaks)):
    f.write(str(peaks_y[i])+" "+str(peaks[i])+" 0 \n")
f.close()



ir=np.loadtxt('VCD.txt',usecols=(0,1,))
idx = (ir[:,0]<2000) & (ir[:,0]>500)
ir=ir[idx]
ir[:,1]=ir[:,1]/np.max(np.abs(ir[:,1]))
idx = np.argsort(ir[:,0])
ir = ir[idx]
#idx = np.argsort(peaks[:,0])
xdat = ir[:,0]
ydat = ir[:,1]
peaks = []
peaks_y = []
for i in range(1,len(ir)-1):
    if(abs(ir[i-1,1])<abs(ir[i,1])>abs(ir[i+1,1])): #and ir[i,0]<2000):
        peaks.append(ir[i,0])
        peaks_y.append(ir[i,1])
f=open('peaks_mod_exp_vcd.txt','w')
for i in range(len(peaks)):
    f.write(str(peaks_y[i])+" "+str(peaks[i])+" 1 \n")
f.close()

ir=np.loadtxt('vcd_theo.txt',usecols=(0,1,))
idx = (ir[:,0]<2000) & (ir[:,0]>500)
ir=ir[idx]
ir[:,1]=ir[:,1]/np.max(np.abs(ir[:,1]))
idx = np.argsort(ir[:,0])
ir = ir[idx]
#idx = np.argsort(peaks[:,0])
xdat = ir[:,0]
ydat = ir[:,1]
peaks = []
peaks_y = []
for i in range(1,len(ir)-1):
    if(abs(ir[i-1,1])<abs(ir[i,1])>abs(ir[i+1,1])): #and ir[i,0]<2000):
        peaks.append(ir[i,0])
        peaks_y.append(ir[i,1])
f=open('peaks_mod_theo_vcd.txt','w')
for i in range(len(peaks)):
    f.write(str(peaks_y[i])+" "+str(peaks[i])+" 1 \n")
f.close()

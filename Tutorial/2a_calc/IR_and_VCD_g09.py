import os
import numpy as np
import matplotlib.pyplot as py
files=os.listdir('.')
freq=[]
inten=[]
vcd=[]
gibbs=[]
for fi in files:
    if(fi.endswith('.log')):
        name_tmp=fi.split('.')[0]
        freq_tmp=[]
        inten_tmp=[]
        vcd_tmp=[]
        f=open(fi,'r')
        for line in f:
            if('Frequencies --' in line):
                tmp=line.split()[2:]
                for i in range(len(tmp)): 
                    freq_tmp.append(float(tmp[i]))
            elif('IR Inten    --' in line):
                tmp=line.split()[3:]
                for i in range(len(tmp)): 
                    inten_tmp.append(float(tmp[i]))
            elif('Rot. str.   --' in line):
                tmp=line.split()[3:]
                for i in range(len(tmp)):
                    vcd_tmp.append(float(tmp[i]))
            elif('Sum of electronic and thermal Free Energies=' in line):
                gibbs_energy=float(line.split()[-1])
        freq_tmp=np.asarray(freq_tmp)
        inten_tmp=np.asarray(inten_tmp)
        vcd_tmp=np.asarray(vcd_tmp)
        if((freq_tmp >= 0).all() and len(freq_tmp)>0):
            freq.append(freq_tmp)
            inten.append(inten_tmp)
            vcd.append(vcd_tmp)
            gibbs.append(gibbs_energy)
freq=np.asarray(freq)
inten=np.asarray(inten)
vcd=np.asarray(vcd)
gibbs=np.asarray(gibbs)
gibbs-=np.min(gibbs)
gibbs*=2625.50
gibbs-=np.min(gibbs)
print(gibbs)
gibbs,idx=np.unique(gibbs,return_index=True)
#idx=np.argsort(gibbs)
#gibbs=gibbs[idx]
print(gibbs)
boltz=np.exp(-gibbs/(0.008314*298.15))
Z=np.sum(boltz)
freq=freq[idx]
vcd=vcd[idx]
inten=inten[idx]
x=np.arange(0,2000)
w_ir=12
w_vcd=12
a=(x[:,np.newaxis,np.newaxis]-freq[np.newaxis,:,:])/(w_ir/2)
ir_spectra=(inten[np.newaxis,:,:]/(1+a**2))
ir_spectra=np.sum(ir_spectra,axis=-1)
y=np.sum(ir_spectra*boltz[np.newaxis,:]/Z,axis=-1)[:,np.newaxis]

a=(x[:,np.newaxis,np.newaxis]-freq[np.newaxis,:,:])/(w_vcd/2)
vcd_spectra=(x[:,np.newaxis,np.newaxis]*vcd[np.newaxis,:,:]/(1+a**2))
vcd_spectra=np.sum(vcd_spectra,axis=-1)
y_vcd=np.sum(vcd_spectra*boltz[np.newaxis,:]/Z,axis=-1)[:,np.newaxis]
print(y_vcd.shape)

ir_spectrum=np.concatenate([x[:,np.newaxis],y],axis=-1)
vcd_spectrum=np.concatenate([x[:,np.newaxis],y_vcd],axis=-1)
np.savetxt('IR_theo.txt',ir_spectrum)
np.savetxt('vcd_theo.txt',vcd_spectrum)

py.plot(ir_spectrum[:,0],ir_spectrum[:,1]/np.max(ir_spectrum[1000:2000,1]))
py.xlim(1000,2000)
py.ylim(0,1)
py.show()


py.plot(vcd_spectrum[:,0],vcd_spectrum[:,1]/np.max(np.abs(vcd_spectrum[1000:2000,1])))
py.xlim(1000,2000)
py.ylim(1,-1)
py.show()





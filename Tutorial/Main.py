import numpy as np
import Algorithm as nd
import matplotlib.pyplot as py
from Settings import Settings
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.interpolate import interp1d
def Lorentzian(freq,inten):
    settings = Settings()
    s,cutoff_absolute,cutoff,sigma_1,sigma_2,use_vcd,x_min,x_max=settings.get()
    x=np.arange(x_min,x_max,1)
    list_append=[]
    for i in range(len(freq)):
        t=((x-freq[i])/(12/2))**2
        L=inten[i]/(1+t)
        list_append.append(L)
    list_append=np.asarray(list_append)
    y=np.sum(list_append,axis=0)
    return x,y



def run():
    settings = Settings()
    s,cutoff_absolute,cutoff,sigma_1,sigma_2,use_vcd,x_min,x_max=settings.get()
    directory=settings.get_directory()
    x_range=np.arange(x_min,x_max)
    for invers in range(-1,2,2):
        vcd=''
        theo_peaks = np.loadtxt(directory+'/peaks_mod_theo.txt',usecols=(0,1,2,))
        theo_peaks[:,1]=theo_peaks[:,1]*s ##Mind that [:,1] is x coordinate
        val = (theo_peaks[:,1]<x_max) & (theo_peaks[:,1]>x_min)
        theo_peaks[:,0] = theo_peaks[:,0]/np.max(theo_peaks[val,0])
        if(use_vcd):
            try:
                theo_peaks_vcd = np.loadtxt(directory+'/peaks_mod_theo_vcd.txt',usecols=(0,1,2,))
                theo_peaks_vcd[:,1]=theo_peaks_vcd[:,1]*s
                theo_peaks_vcd[:,0]=theo_peaks_vcd[:,0]*invers
                val = (theo_peaks_vcd[:,1]<x_max) & (theo_peaks_vcd[:,1]>x_min)
                theo_peaks_vcd[:,0]=theo_peaks_vcd[:,0]/np.max(np.abs(theo_peaks_vcd[val,0]))
                theo_peaks=np.concatenate([theo_peaks,theo_peaks_vcd],axis=0)
                idx=np.argsort(theo_peaks[:,1])

                theo_peaks=theo_peaks[idx]
                vcd='_vcd_'
                print("VCD SPECTRUM AVAILABLE")
            except:
                print("NO VCD SPECTRUM AVAILABLE")
                pass

        exp_peaks = np.loadtxt(directory+'/peaks_mod_exp.txt',usecols=(0,1,2,))
        val = (exp_peaks[:,1]<x_max) & (exp_peaks[:,1]>x_min)
        exp_peaks[:,0]=exp_peaks[:,0]/np.max(exp_peaks[val,0])
        if(use_vcd):
            try:
                exp_peaks_vcd = np.loadtxt(directory+'/peaks_mod_exp_vcd.txt',usecols=(0,1,2,))
                val = (exp_peaks_vcd[:,1]<x_max) & (exp_peaks_vcd[:,1]>x_min)
                exp_peaks_vcd[:,0]=exp_peaks_vcd[:,0]/np.max(np.abs(exp_peaks_vcd[val,0]))
                exp_peaks=np.concatenate([exp_peaks,exp_peaks_vcd],axis=0)
                idx=np.argsort(exp_peaks[:,1])
                exp_peaks=exp_peaks[idx]
            except:
                print("NO EXPERIMENTAL VCD SPECTRUM AVAILABLE")
                pass
        print("CALL ALIGNMENT")
        Algorithm = nd.Algorithm(theo_peaks,exp_peaks,cutoff=cutoff)
        print("START ALIGNMENT")
        returnvalue,old_freq,freq,inten,vcd_ir_array=Algorithm.Needleman_IR()
        print("ALIGNMENT SUCCEESFUL")
        print('returnvalue',returnvalue)
        print("READ IN EXPERIMENT")
        ir_exp=np.loadtxt(directory+"/IR.txt",usecols=(0,1,))

        f_ir_exp=interp1d(ir_exp[:,0],ir_exp[:,1],kind='cubic')
        val=(ir_exp[:,0]>x_min) & (ir_exp[:,0]<x_max)
        ir_exp[:,1]=ir_exp[:,1]/np.max(ir_exp[val,1])
        print("READ IN THEORETICAL SPECTRUM")
        theo_unshifted=np.loadtxt(directory+'/IR_theo.txt',usecols=(0,1,))

        vcd_ir_array=np.asarray(vcd_ir_array,dtype=int)
        x,y=Lorentzian(freq[vcd_ir_array==0],inten[vcd_ir_array==0])

        val = (theo_unshifted[:,0]*s>x_min) & (theo_unshifted[:,0]*s<x_max)
        py.plot(ir_exp[:,0], ir_exp[:,1],label="EXPERIMENTAL")
        theo_unshifted[:,1]=theo_unshifted[:,1]/np.max(theo_unshifted[val,1])
        py.plot(theo_unshifted[:,0],theo_unshifted[:,1],label="UNSHIFTED")
        
        val = (x<x_max) & (x>x_min)
        y=y/np.max(y[val])
        f_ir_theo=interp1d(x,y,kind='cubic')
        py.plot(x,y,label="SHIFTED")
        py.legend()
        py.xlim(x_max,x_min)
        fig = py.figure(figsize=(4.2,3.6))
        ax = fig.add_subplot(111)
        ax.plot(ir_exp[:,0],ir_exp[:,1],color='black',linewidth=1,label='EXP')
        ax.plot(theo_unshifted[:,0]*s,theo_unshifted[:,1],color='orange',linewidth=1,label='UNSHIFTED')
        ax.plot(x,y,color='red',linewidth=1,label='SHIFTED')
        ax.set_xlim(x_max,x_min)
        ax.set_ylim(-0.1,1.01)
        ax.set_xlabel("$\\tilde{\\nu}$ / cm$^{-1}$",labelpad=0)
        ax.set_ylabel("$I$ / a.u.",labelpad=0)
        #py.show()
        fig.savefig(directory+"_0_"+str(s)+"_"+str(cutoff)+vcd+"_"+str(invers)+".png",dpi=500)
        if(use_vcd):
            try:
                vcd_exp=np.loadtxt(directory+'/VCD.txt')
                val=(vcd_exp[:,0]>x_min) & (vcd_exp[:,0]<x_max)
                vcd_exp[:,1]=vcd_exp[:,1]/np.max(np.abs(vcd_exp[val,1]))

                vcd_theo=np.loadtxt(directory+'/vcd_theo.txt')
                val=(vcd_theo[:,0]>x_min) & (vcd_theo[:,0]<x_max)
                vcd_theo[:,1] = vcd_theo[:,1]/np.max(np.abs(vcd_theo[val,1]))
                x,y=Lorentzian(freq[vcd_ir_array==1],inten[vcd_ir_array==1])
                fig = py.figure(figsize=(4.2,3.6))
                ax = fig.add_subplot(111)
                ax.plot(vcd_exp[:,0],vcd_exp[:,1],color='black',linewidth=1,label='EXP')
                ax.plot(vcd_theo[:,0]*s,invers*vcd_theo[:,1],color='orange',linewidth=1,label='UNSHIFTED')
                val=(x>x_min) & (x<x_max)
                y=y/np.max(np.abs(y[val]))
                f_vcd_theo=interp1d(x,y,kind='cubic')
                f_vcd_exp=interp1d(vcd_exp[:,0],vcd_exp[:,1],kind='cubic')
                ax.plot(x,y,color='red',linewidth=1,label='SHIFTED')
                
                ax.set_xlim(x_max,x_min)
                ax.set_ylim(-1.01,1.01)
                ax.set_xlabel("$\\tilde{\\nu}$ / cm$^{-1}$",labelpad=0)
                ax.set_ylabel("$I$ / a.u.",labelpad=0)
                fig.savefig(directory+"_1_"+str(s)+"_"+str(cutoff)+vcd+"_"+str(invers)+".png",dpi=500)
            except:
                print("VCD NOT PLOTTED SINCE NOT EXISTENT")
                pass
        ##COMPUTE METRICS
        print(invers,"Alignment Score", -returnvalue)
        print(invers,"IR pearson number",pearsonr(f_ir_theo(x_range),f_ir_exp(x_range))[0])
        print(invers,"IR spearman number",spearmanr(f_ir_theo(x_range),f_ir_exp(x_range))[0])
        try:
            print(invers,"VCD pearson number",pearsonr(f_vcd_theo(x_range),f_vcd_exp(x_range))[0])
            print(invers,"VCD spearman number",spearmanr(f_vcd_theo(x_range),f_vcd_exp(x_range))[0])
        except:
            print("VCD metrics cannot be computed")
            pass
    return 0


EXIT_CODE = run()
if(EXIT_CODE == 0):
    print("sucess")
    exit()
else:
    print("fail")
    exit()

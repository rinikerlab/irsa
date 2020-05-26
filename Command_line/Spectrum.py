import numpy as np
import sys
import os
import matplotlib.pyplot as py
from scipy.stats import pearsonr
""" This file implements the class Spectrum:
    class Spectrum contains all information about the unshifted, boltzmann weighted and normalized (at T=300 K) theoretical
    spectrum and experimental spectrum, and
    about the shifted, boltzmann weighted theoretical spectrum.
    It saves the theoretical spectrum to a png file. """

class Spectrum:
    def __init__(self,exp,freq,E,T,lower_bound,higher_bound,w1,w2,resolution_exp,Dipol,out):
        """init function. Sets the theoretical and experimental spectra, and compute all necessary information """
        self.out = out
        self.resolution_exp = resolution_exp
        self.exp = exp
        self.freq = freq
        self.E = E
        self.T = T
        self.u = lower_bound #lower bound for the match considered
        self.h = higher_bound #higher bound for the match considered
        self.w1 = w1 #Lorentzbandwidth unshifted
        self.w2 = w2 #Lorentzbandwidth shifted
        self.Dipol = Dipol #Whether the pickle file was saved as DIPOL or IR INTENSITIES
        self.compute_Z() #compute Boltzmannfactor
        self.Boltzmann_weight() #compute weighted Spectrum
        self.normalize() #normalize the boltzmann weighted spectrum
        self.normalize_exp() #normalize the experimental spectrum
        self.theo_peaks()
        self.exp_peaks()
        print("Object spectrum successfully created")
    def compute_Z(self): 
        """Compute Z in N_i/N = \sum_i e^(E_i/(RT))/Z"""
        Z_array_tmp = np.exp(-self.E/(self.T*8.314*10**(-3)))
        self.Z = np.sum(Z_array_tmp,axis=-1)
    def Boltzmann_weight(self):
        """Computes the Boltzmann weighted spectrum"""
        self.spectrum = np.zeros((len(self.freq),2000,2)) ## contains the broadened spectrum of each conformer
        self.spectrum_boltzmann = np.zeros((2000,2))    ## contains the boltzmannweighted spectrum
        ## x_axis of theoretical spectrum
        x_value = np.arange(0,2000,1)
        self.spectrum[:,:,0] = x_value
        ##Lorentz broadening for all conformers
        x = (self.freq[:,:,np.newaxis,0]-x_value[np.newaxis,:])/(self.w1/2) ## first factor for Lorentz broadening
        if(self.Dipol == 0):
            self.spectrum[:,:,1] = np.sum((self.freq[:,:,1,np.newaxis]/(1+x**2)[:,:])*x_value/(2.236*10**(-39)),axis=1)
        else: ## IR INTENSITIES
            self.spectrum[:,:,1] = np.sum((self.freq[:,:,1,np.newaxis]/(1+x**2)[:,:]),axis=1)
        ##Create the boltzmann weighted spectrum
        self.spectrum_boltzmann[:,0] = x_value ##x_axis 0..2000
        if(len(self.E)>1):
            self.spectrum_boltzmann[:,1] = np.sum(self.spectrum[:,:,1]*self.E[:,np.newaxis]/self.Z,axis=0)
        else:
            self.spectrum_boltzmann[:,1] = self.spectrum[0,:,1]#*self.E[0]/self.Z
    def normalize(self):
        """Normalizes the theoretical spectrum"""
        self.spectrum_boltzmann[:,1] = self.spectrum_boltzmann[:,1]/np.max(self.spectrum_boltzmann[self.u:self.h,1])
    def normalize_exp(self):
        """Normalizes the experimental spectrum. """
        u_exp = 0
        h_exp = 0
        for i in range(len(self.exp)):
            if(abs(self.exp[i,0]-self.u) < self.resolution_exp and u_exp == 0 ):
                u_exp = i
            if(abs(self.exp[i,0]-self.h) < self.resolution_exp and h_exp == 0 ):
                h_exp = i
        self.exp[:,1]=self.exp[:,1]/np.max(self.exp[u_exp:h_exp,1])
    def theo_peaks(self):
        """Automatically pick peaks in the theoretical spectrum. Simplest criterium one could think off. Place for
           improvement/manually selection if it fails (NOT YET IMPLEMENTED IN THIS VERSION, REQUEST) """
        theo_x = []
        theo_y = []
        for i in range(1,len(self.spectrum_boltzmann)-1):
            if(self.spectrum_boltzmann[i,0]>self.u and self.spectrum_boltzmann[i,0]<self.h):
                if(self.spectrum_boltzmann[i-1,1]<self.spectrum_boltzmann[i,1]>self.spectrum_boltzmann[i+1,1]):
                    theo_x.append(self.spectrum_boltzmann[i,0])
                    theo_y.append(self.spectrum_boltzmann[i,1])
        self.theo_peaks = np.zeros((len(theo_x),2))
        self.theo_peaks[:,0] = np.asarray(theo_x)
        self.theo_peaks[:,1] = np.asarray(theo_y)
    def exp_peaks(self):
        """Automatically pick peaks in the experimental spectrum. Simplest criterium one could think off. Place for
           improvement/manually selection if it fails (NOT YET IMPLEMENTED IN THIS VERSION, REQUEST)"""
        exp_x = []
        exp_y = []
        for i in range(1,len(self.exp)-1):
            if(self.exp[i,0]>self.u and self.exp[i,0]<self.h):
                if(self.exp[i-1,1]<self.exp[i,1]>self.exp[i+1,1]):
                    exp_x.append(self.exp[i,0])
                    exp_y.append(self.exp[i,1])
        self.exp_peaks = np.zeros((len(exp_x),2))
        self.exp_peaks[:,0] = np.asarray(exp_x)
        self.exp_peaks[:,1] = np.asarray(exp_y)
    def set_new_freq(self,freq_old,freq_new,inten_new):
        self.freq_new = freq_new
        self.inten_new = inten_new
        self.freq_old = freq_old # THIS LINE IS REDUNDANT
    """PLOT FUNCTION OF THE UNSHFITED PLOT WITH THE ASSIGNMENTS MADE"""
    def plot_unshifted(self):
        self.spectrum_boltzmann[:,1] = self.spectrum_boltzmann[:,1]/np.max(self.spectrum_boltzmann[self.u:self.h,1])
        py.plot(self.spectrum_boltzmann[:,0],self.spectrum_boltzmann[:,1],color="red")
        py.plot(self.exp[:,0],self.exp[:,1],color="black")
        py.plot()
        for i in range(len(self.freq_old)):
            py.plot([self.freq_old[i],self.freq_new[i]],[self.inten_new[i],self.inten_new[i]],color="blue")
        py.plot(self.theo_peaks[:,0],self.theo_peaks[:,1],".",color="black")
        py.plot(self.exp_peaks[:,0],self.exp_peaks[:,1],"x",color="black")
        py.xlim(self.h,self.u)
        py.ylim(0,1.02)
        py.savefig(self.out+"unshifted.png",dpi=400)
        py.cla()
        py.clf()
        py.close()
    def create_shifted(self):
        """Computes the Boltzmann weighted spectrum"""
        self.freq_shifted = np.zeros((1,len(self.freq_new),2))
        self.freq_shifted[0,:,0] = np.asarray(self.freq_new)
        self.freq_shifted[0,:,1] = np.asarray(self.inten_new)
        self.spectrum_shifted = np.zeros((1,2000,2)) ## contains the broadened spectrum of each conformer
        self.spectrum_bc = np.zeros((2000,2))
        ## x_axis of theoretical spectrum
        x_value = np.arange(0,2000,1)
        self.spectrum_shifted[:,:,0] = x_value
        ##Lorentz broadening for all conformers
        x = (self.freq_shifted[:,:,np.newaxis,0]-x_value[np.newaxis,:])/(self.w2/2) ## first factor for Lorentz broadening
        self.spectrum_shifted[:,:,1] = np.sum((self.freq_shifted[:,:,1,np.newaxis]/(1+x**2)[:,:]),axis=1)
        ##Create the boltzmann weighted spectrum
        self.spectrum_bc[:,0] = self.spectrum_shifted[0,:,0] ##x_axis 0..2000
        self.spectrum_bc[:,1] = self.spectrum_shifted[0,:,1]
        """IMPLEMENTATION OF ALL THE GETTER FUNCTIONS"""
    def get_lower_bound(self):
        return self.u
    def get_higher_bound(self):
        return self.h
    def get_theo_peaks(self):
        return self.theo_peaks
    def get_exp_peaks(self):
        return self.exp_peaks 
    def plot(self):
        """PLOT ALIGNED SPECTRUM"""
        self.spectrum_bc[:,1] = self.spectrum_bc[:,1]/np.max(self.spectrum_bc[self.u:self.h,1])
        py.plot(self.spectrum_bc[:,0],self.spectrum_bc[:,1],color="red")
        py.plot(self.exp[:,0],self.exp[:,1],color="black")
        py.xlim(self.h,self.u)
        py.ylim(0,1.02)
        py.savefig(self.out+"aligned.png",dpi=800)
        py.cla()
        py.clf()
        py.close()
    def integrate(self,number_of_bins=150):
        """INTEGRATE THE EXPERIMENTAL AND THEORETICAL SPECTRUM AND COMPUTES THE PEARSON COEFFICIENT (OTHER CORRELATIONS
        COEFFICIENTS MIGHT BE USED AS WELL) """
        spaceing = (int)((self.h-self.u)/number_of_bins)
        integ_theo = np.zeros((number_of_bins,2))
        integ_exp = np.zeros((number_of_bins,2))
        for i in range(number_of_bins):
            integ_theo[i,0]=spaceing*i+self.u
            counter = 0
            for j in range(len(self.spectrum_bc)):
                if(spaceing*i+self.u<=self.spectrum_bc[j,0]<spaceing*(i+1)+self.u):
                    integ_theo[i,1]+=self.spectrum_bc[j,1]
                    counter+=1
            if(counter>0):
                integ_theo[i,1]/=counter
        for i in range(number_of_bins):
            integ_exp[i,0]=spaceing*i+self.u
            counter = 0
            for j in range(len(self.exp)):
                if(spaceing*i+self.u<=self.exp[j,0]<spaceing*(i+1)+self.u):
                    integ_exp[i,1]+=self.exp[j,1]
                    counter+=1
            if(counter>0):
                integ_exp[i,1]/=counter
        p = pearsonr(integ_exp[:,1],integ_theo[:,1])[0]
        return p

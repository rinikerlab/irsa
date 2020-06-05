import numpy as np
import sys
import os
import matplotlib.pyplot as py
from scipy.stats import pearsonr
import global_ as g
""" This file implements the class Spectrum:
    class Spectrum contains all information about the unshifted, boltzmann weighted and normalized (at T=300 K) theoretical
    spectrum and experimental spectrum, and
    about the shifted, boltzmann weighted theoretical spectrum.
    It saves the theoretical spectrum to a png file. """

class Spectrum:
    def __init__(self):
        """init function. Sets the theoretical and experimental spectra, and compute all necessary information """
        print("Object spectrum successfully created")
    def compute_Z(self): 
        """Compute Z in N_i/N = \sum_i e^(E_i/(RT))/Z"""
        g.E -= np.min(g.E)
        g.Z_array = np.exp(-g.E/(g.values['T']*8.314*10**(-3)))
        g.Z = np.sum(g.Z_array,axis=-1)
        g.Z_array_c = g.Z_array/g.Z
    def Boltzmann_weight_IR(self):
        self.compute_Z()
        """Computes the Boltzmann weighted spectrum"""
        g.spectrum = np.zeros((len(g.theo_ir),2000,2)) ## contains the broadened spectrum of each conformer
        g.spectrum_boltzmann = np.zeros((2000,2))    ## contains the boltzmannweighted spectrum
        ## x_axis of theoretical spectrum
        x_value = np.arange(0,2000,1)
        g.spectrum[:,:,0] = x_value
        ##Lorentz broadening for all conformers
        x = (g.theo_ir[:,:,np.newaxis,0]-x_value[np.newaxis,:])/(g.values['w']/2) ## first factor for Lorentz broadening
        if(g.Dipol == 0):
            g.spectrum[:,:,1] = np.sum((g.theo_ir[:,:,1,np.newaxis]/(1+x**2)[:,:])*x_value/(2.236*10**(-39)),axis=1)
        else: ## IR INTENSITIES
            g.spectrum[:,:,1] = np.sum((g.theo_ir[:,:,1,np.newaxis]/(1+x**2)[:,:]),axis=1)
        ##Create the boltzmann weighted spectrum
        g.spectrum_boltzmann[:,0] = x_value ##x_axis 0..2000
        if(len(g.E)>1):
            g.spectrum_boltzmann[:,1] = np.sum(g.spectrum[:,:,1]*g.Z_array_c[:,np.newaxis],axis=0)
        else:
            g.spectrum_boltzmann[:,1] = g.spectrum[0,:,1]
        g.spectrum_boltzmann[:,1] = g.spectrum_boltzmann[:,1]/np.max(np.abs(g.spectrum_boltzmann[g.values['lb']:g.values['hb'],1]))
    def Boltzmann_weight_VCD(self):
        self.compute_Z()
        """Computes the Boltzmann weighted spectrum"""
        g.spectrum_vcd = np.zeros((len(g.theo_ir),2000,2)) ## contains the broadened spectrum of each conformer
        g.spectrum_boltzmann_vcd = np.zeros((2000,2))    ## contains the boltzmannweighted spectrum
        ## x_axis of theoretical spectrum
        x_value = np.arange(0,2000,1)
        g.spectrum_vcd[:,:,0] = x_value
        ##Lorentz broadening for all conformers
        x = (g.theo_vcd[:,:,np.newaxis,0]-x_value[np.newaxis,:])/(g.values['w']/2) ## first factor for Lorentz broadening
        g.spectrum_vcd[:,:,1] = np.sum((g.theo_vcd[:,:,1,np.newaxis]/(1+x**2)[:,:])*x_value/(2.236*10**(-39)),axis=1) ## VCD ROT
        ##Create the boltzmann weighted spectrum
        g.spectrum_boltzmann_vcd[:,0] = x_value ##x_axis 0..2000
        if(len(g.E)>1):
            g.spectrum_boltzmann_vcd[:,1] = np.sum(g.spectrum_vcd[:,:,1]*g.Z_array_c[:,np.newaxis],axis=0)
        else:
            g.spectrum_boltzmann_vcd[:,1] = g.spectrum_vcd[0,:,1]
        g.spectrum_boltzmann_vcd[:,1] = g.spectrum_boltzmann_vcd[:,1]/np.max(np.abs(g.spectrum_boltzmann_vcd[g.values['lb']:g.values['hb'],1]))
    def VCD_shifted(self):
        x_value = np.arange(0,2000,1)
        g.VCD_shifted = np.zeros((2000,2))
        g.VCD_shifted[:,0] = x_value
        x = (np.asarray(np.asarray(g.freq_new)[np.newaxis,:])-x_value[:,np.newaxis])/(g.values['w']/2) ## first factor for Lorentz broadening
        g.VCD_shifted[:,1] = np.sum((np.asarray(g.inten_VCD_new)[np.newaxis,:]/(1+x**2)),axis=1)
    def IR_shifted(self):
        x_value = np.arange(0,2000,1)
        g.IR_shifted = np.zeros((2000,2))
        g.IR_shifted[:,0] = x_value ## (2000,30)
        x = (np.asarray(np.asarray(g.freq_new)[np.newaxis,:])-x_value[:,np.newaxis])/(g.values['w']/2) ## first factor for Lorentz broadening ((2000,30))
        g.IR_shifted[:,1] = np.sum(np.asarray(g.inten_new)[np.newaxis,:]/(1+x**2),axis=1)
    def integrate(self,number_of_bins=150):
        def i(data_set):
            spaceing = (int)((g.values['hb']-g.values['lb'])/number_of_bins)
            integ = np.zeros((number_of_bins,2))
            for i in range(number_of_bins):
                integ[i,0]=spaceing*i+g.values['lb']
                counter = 0
                for j in range(len(data_set)):
                    if(spaceing*i+g.values['lb']<=data_set[j,0]<spaceing*(i+1)+g.values['lb']):
                        integ[i,1]+=data_set[j,1]
                        counter=counter+1
                if(counter>0):
                    integ[i,1]=integ[i,1]/counter
            return integ
        IR_theo_integ = i(g.IR_shifted)
        IR_exp_integ = i(g.exp_ir)
        p_ir = pearsonr(IR_theo_integ[:,1],IR_exp_integ[:,1])[0]
        try:
            VCD_theo_integ = i(g.VCD_shifted)
            VCD_exp_integ = i(g.exp_vcd)
            p_vcd = pearsonr(VCD_theo_integ[:,1],VCD_exp_integ[:,1])[0]
            return p_ir,p_vcd
        except: 
            pass
        return p_ir


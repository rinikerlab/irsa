import numpy as np
from Settings import Settings
import os

"""THIS CLASS IMPLEMENTS THE SHIFTING ALGORITHM. CARE ABOUT THE VALUES IN THE INIT FUNCTION. DEPEND ON THE ELECTRONIC
   STRUCTURE THEORY CHOSEN"""
class Algorithm:
    def __init__(self,theo_peaks,exp_peaks,cutoff=0.01):
        """SPECTRUM INFORMATION"""
        print("Initialization ND")
        settings = Settings()
        self.s,self.cutoff_absolute,self.cutoff,self.sigma_1,self.sigma_2,self.use_vcd,self.u,self.h=settings.get()
        self.cutoff=cutoff
        self.theo_peaks=theo_peaks
        self.exp_peaks=exp_peaks

        print("Initialization SUCCESSFUL")
    def Diagonal_IR(self,freq_i,inten_i,exp_freq_j,exp_inten_j,bond_l,bond_h,exp_vcd=0,inten_vcd=0):
        """COMPUTE THE SCORES FOR EACH PEAK COMBINATION DYNAMICALLY"""
        value=0.01
        if(inten_vcd==exp_vcd):
            if(((min(abs(1-exp_freq_j/freq_i),abs(1-freq_i/exp_freq_j))<self.cutoff and exp_freq_j > self.u and exp_freq_j < self.h and self.cutoff_absolute==False) or (self.cutoff_absolute==True and abs(exp_freq_j-freq_i)<self.cutoff)) and bond_l< exp_freq_j < bond_h and bond_l < freq_i < bond_h):
                inten_contrib = np.exp(-0.5*(exp_freq_j/freq_i)**2/self.sigma_1)
                freq_contrib = np.exp(-0.5*(abs(1-min(inten_i/exp_inten_j,exp_inten_j/inten_i)))**2/self.sigma_2)
                value = (-1)*freq_contrib*inten_contrib#-inten_contrib#*eta_contrib
                if(inten_vcd==1):
                    value*=0.1
        return value
    def Backtrace_IR(self,p_mat,al_mat,n,m,freq_i,inten_i,exp_freq_j,exp_inten_i,bond_l,bond_h,vcd): #n theoretical, m experimental
        #freq, inten, old_freq, new_sigma = self.Backtrace_IR(p_mat,al_mat,n,m,freq,inten,exp_freq,sigma,bond_l,bond_h,exp_sigma,vcd=vcd)
        """BACKTRACE THE NEEDLEMAN ALGORITHM"""
        #ir_unshifted_x, ir_unshifted_y, vcd_unshifted_x, vcd_unshifted_y = self.Spectrum.get_non_matched()
        new_freq = []
        new_freq_VCD = [] ## only VCD peaks
        old_freq = []
        new_inten = []
        new_inten_vcd = []
        non_matched = []
        matched_freq = []
        vcd_ir_array = []
        n=n-1
        m=m-1
        current_scaling_factor = 1
        factors = []
##FIRST, SHIFT ALL THE PEAKS WHICH ARE MATCHED BY THE ALGORITHM
        while(True):
            if(p_mat[n,m]=="D"):
                new_freq.append(exp_freq_j[m-1])
                old_freq.append(freq_i[n-1])
                new_inten.append(inten_i[n-1])
                vcd_ir_array.append(vcd[n-1])
                current_scaling_factor = exp_freq_j[m-1]/freq_i[n-1]
                matched_freq.append(n-1)
                factors.append(current_scaling_factor)
                n=n-1
                m=m-1
            elif(p_mat[n,m]=="V"):
                non_matched.append(n-1)
                n=n-1
            elif(p_mat[n,m]=="H"):
                m=m-1
            else:
                break
##SECOND, MATCH ALL PEAKS, WHICH ARE NOT MATCHED
        for i in range(len(non_matched)):
            closest_distance = 9999
            matched_to = 0
            sf = 1
            for j in range(len(matched_freq)):
                dis=abs(freq_i[non_matched[i]]-freq_i[matched_freq[j]])
                if(dis<closest_distance):
                    closest_distance = dis
                    sf = factors[j]
            new_freq.append(freq_i[non_matched[i]]*sf)
            vcd_ir_array.append(vcd[non_matched[i]])
            old_freq.append(freq_i[non_matched[i]])
            new_inten.append(inten_i[non_matched[i]])
        return np.asarray(new_freq),np.asarray(new_inten),np.asarray(old_freq),np.asarray(vcd_ir_array)
    def Pointer(self,di,ho,ve):
        """POINTER TO CELL IN THE TABLE"""
        pointer = min(di,min(ho,ve))
        if(di==pointer):
            return "D"
        elif(ho==pointer):
            return "H"
        else:
            return "V"
    def Needleman_IR(self):
        """NEEDLEMAN ALGORITHM FOR IR"""
        freq = self.theo_peaks[:,1]
        inten = self.theo_peaks[:,0]
        vcd = self.theo_peaks[:,2]
        exp_freq = self.exp_peaks[:,1]
        exp_inten = self.exp_peaks[:,0]
        exp_inten_vcd = self.exp_peaks[:,2]
        bond_l = self.u
        bond_h = self.h
        n = len(freq)+1
        m = len(exp_freq)+1
        norm = 1
        al_mat = np.zeros((n,m))
        p_mat = np.zeros((n,m),dtype='U25') #string
        for i in range(1,n):
            al_mat[i,0] = al_mat[i-1,0]#+0.01#self.dummy_0 # BOUND SOLUTION, VALUE MIGHT BE CHANGED
            p_mat[i,0] = 'V'
        for i in range(1,m):
            al_mat[0,i] = al_mat[0,i-1]#+0.01##+self.dummy_1
            p_mat[0,i] = 'H'
        p_mat[0,0]="S"
        normalize = 0
        for i in range(1,n): #theoretical
            for j in range(1,m): #experimental
                di = self.Diagonal_IR(freq[i-1],inten[i-1],exp_freq[j-1],exp_inten[j-1],bond_l=self.u,bond_h=self.h,exp_vcd=exp_inten_vcd[j-1],inten_vcd=vcd[i-1])
                di = al_mat[i-1,j-1]+di
                ho = al_mat[i,j-1]
                ve = al_mat[i-1,j]
                al_mat[i,j] = min(di,min(ho,ve))
                p_mat[i,j] = self.Pointer(di,ho,ve)
        freq, inten, old_freq, vcd = self.Backtrace_IR(p_mat,al_mat,n,m,freq,inten,exp_freq,exp_inten,bond_l=self.u,bond_h=self.h,vcd=vcd)
        returnvalue = al_mat[n-1,m-1]#/(n+m) ##ORIGINALLY WE DIVIDED BY THE NUMBER OF THEORETICAL PEAKS
                                           ##HOWEVER, WE FOUND THIS TOO INCONVIENT, SINCE IT MAKES THE DEPENDENCE ON THE
                                           ##PURE NUMBERS TOO LARGE
        return returnvalue, old_freq, freq, inten,np.asarray(vcd)


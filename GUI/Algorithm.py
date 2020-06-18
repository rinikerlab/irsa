import numpy as np
import Spectrum
import os
import global_ as g

"""THIS CLASS IMPLEMENTS THE SHIFTING ALGORITHM. CARE ABOUT THE VALUES IN THE INIT FUNCTION. DEPEND ON THE ELECTRONIC
   STRUCTURE THEORY CHOSEN"""
class Algorithm:
    def __init__(self):
        """SPECTRUM INFORMATION"""
        print("ALGORITHM CREATED")
        self.Spectrum = g.Spectrum
        self.u = g.values["lb"]
        self.h = g.values["hb"]
        self.exp_peaks = g.exp_peaks 
        self.theo_peaks = g.theo_peaks
        self.n = len(self.theo_peaks)# number theo peaks
        self.m = len(self.exp_peaks)# number exp peaks
        """ALGORITHM SETTINGS"""
        self.sigma_0 = g.values['s0']
        self.sigma_1 = g.values['s1']
        self.cutoff = g.values['c']
        self.dummy_0 = 0
        self.dummy_1 = 0
        self.mu = g.values['mu']
        if(g.set_VCD==True):
            self.theo_peaks_vcd = np.zeros((len(g.theo_peaks_x),2))
            self.theo_peaks_vcd[:,0] = np.asarray(g.theo_peaks_x)
            self.theo_peaks_vcd[:,1] = np.asarray(g.peak_list_VCD_y_theo)
        print("ASSIGNED 2")
    def Diagonal_IR(self,freq_i,inten_i,exp_freq_j,exp_inten_j,bond_l,bond_h,n,m):
        """COMPUTE THE SCORES FOR EACH PEAK COMBINATION DYNAMICALLY"""
        if(abs(exp_freq_j-freq_i*self.mu)<self.cutoff):
            value = (-1)*np.exp(-0.5*(min(inten_i/exp_inten_j,exp_inten_j/inten_i)-1)**2/self.sigma_0**2)*np.exp(-0.5*(exp_freq_j/(freq_i)-self.mu)**2/self.sigma_1**2)
        else:
            value = 1
        return value
    def Diagonal_VCD(self,freq_i,inten_i,inten_i_VCD,exp_freq_j,exp_inten_j,exp_inten_j_VCD,bond_l,bond_h,n,m):
        if(abs(exp_freq_j-freq_i*self.mu)<self.cutoff):
            value = (-1)*exp_inten_j_VCD*inten_i_VCD*np.exp(-0.5*(min(inten_i/exp_inten_j,exp_inten_j/inten_i)-1)**2/self.sigma_0**2)*np.exp(-0.5*(exp_freq_j/(freq_i)-self.mu)**2/self.sigma_1**2)
        else:
            value = 1
        return value
    def Backtrace_IR(self,p_mat,al_mat,n,m,freq_i,inten_i,exp_freq_j,bond_l,bond_h): #n theoretical, m experimental
        """BACKTRACE THE NEEDLEMAN ALGORITHM"""
        new_freq = []
        old_freq = []
        new_inten = []
        non_matched_freq = []
        matched_freq = []
        non_matched_inten = []
        n=n-1
        m=m-1
        current_scaling_factor = 1
        factors = []
        while(True):
            if(p_mat[n,m]=="D"):
                new_freq.append(exp_freq_j[m-1])
                old_freq.append(freq_i[n-1])
                new_inten.append(inten_i[n-1])
                current_scaling_factor = exp_freq_j[m-1]/freq_i[n-1]
                matched_freq.append(n-1)
                factors.append(current_scaling_factor)
                n=n-1
                m=m-1
            elif(p_mat[n,m]=="V"):
                non_matched_inten.append(n-1)
                non_matched_freq.append(n-1)
                n=n-1
            elif(p_mat[n,m]=="H"):
                m=m-1
            else:
                break
        for i in range(len(non_matched_freq)):
            closest_distance = 9999
            matched_to = 0
            sf = 1
            for j in range(len(matched_freq)):
                dis=abs(freq_i[non_matched_freq[i]]-freq_i[matched_freq[j]])
                if(dis<closest_distance):
                    closest_distance = dis
                    sf = factors[j]
            new_freq.append(freq_i[non_matched_freq[i]]*sf)
            old_freq.append(freq_i[non_matched_freq[i]])
            new_inten.append(inten_i[non_matched_freq[i]])
        return new_freq,new_inten,old_freq
    def Backtrace_VCD(self,p_mat,al_mat,n,m,freq_i,inten_i,inten_i_VCD,exp_freq_j,bond_l,bond_h): #n theoretical, m experimental
        new_freq = []
        old_freq = []
        new_inten = []
        new_inten_VCD = []
        non_matched_freq = []
        matched_freq = []
        non_matched_inten = []
        non_matched_inten_VCD = []
        n=n-1
        m=m-1
        current_scaling_factor = 1
        factors = []
        while(True):
            if(p_mat[n,m]=="D"):
                old_freq.append(freq_i[n-1])
                new_freq.append(exp_freq_j[m-1])    
                new_inten.append(inten_i[n-1])
                new_inten_VCD.append(inten_i_VCD[n-1])
                current_scaling_factor = exp_freq_j[m-1]/freq_i[n-1]
                matched_freq.append(n-1)
                factors.append(current_scaling_factor)
                n=n-1
                m=m-1
            elif(p_mat[n,m]=="V"):
                non_matched_inten.append(n-1)
                non_matched_inten_VCD.append(n-1)
                non_matched_freq.append(n-1)
                n=n-1
            elif(p_mat[n,m]=="H"):
                m=m-1
            else:
                break
        for i in range(len(non_matched_freq)):
            closest_distance = 9999
            matched_to = 0
            sf = 1
            for j in range(len(matched_freq)):
                dis=abs(freq_i[non_matched_freq[i]]-freq_i[matched_freq[j]])
                if(dis<closest_distance):
                    closest_distance = dis
                    sf = factors[j]
            new_freq.append(freq_i[non_matched_freq[i]]*sf)
            old_freq.append(freq_i[non_matched_freq[i]])
            new_inten.append(inten_i[non_matched_freq[i]])
            new_inten_VCD.append(inten_i_VCD[non_matched_freq[i]])
        return new_freq,new_inten,new_inten_VCD,old_freq
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
        freq = self.theo_peaks[:,0]
        inten = self.theo_peaks[:,1]
        exp_freq = self.exp_peaks[:,0]
        exp_inten = self.exp_peaks[:,1]
        if(g.set_VCD == True):
            inten_VCD = np.asarray(self.theo_peaks_vcd[:,1])
            exp_inten_VCD = np.asarray(g.peak_list_VCD_y)
        bond_l = self.u
        bond_h = self.h
        n = self.n+1
        m = self.m+1
        norm = 1
        al_mat = np.zeros((n,m))
        p_mat = np.zeros((n,m),dtype='U25') #string
        for i in range(1,n):
            al_mat[i,0] = al_mat[i-1,0]+self.dummy_0 # BOUND SOLUTION, VALUE MIGHT BE CHANGED
            p_mat[i,0] = 'V'
        for i in range(1,m):
            al_mat[0,i] = al_mat[0,i-1]+self.dummy_1
            p_mat[0,i] = 'H'
        p_mat[0,0]="S"
        normalize = 0
        for i in range(1,n): #theoretical
            for j in range(1,m): #experimental
                print("hoho")
                if(g.set_VCD == False):
                    di = self.Diagonal_IR(freq[i-1],inten[i-1],exp_freq[j-1],exp_inten[j-1],bond_l,bond_h,n,m)
                else:
    #def Diagonal_VCD(freq_i,inten_i,inten_i_VCD,exp_freq_j,exp_inten_j,exp_inten_j_VCD,bond_l,bond_h,n,m):
    #def Diagonal_VCD(self,freq_i,inten_i,inten_i_VCD,exp_freq_j,exp_inten_j,exp_inten_j_VCD,bond_l,bond_h,n,m):
                    di = self.Diagonal_VCD(freq[i-1],inten[i-1],inten_VCD[i-1],exp_freq[j-1],exp_inten[j-1],exp_inten_VCD[j-1],bond_l,bond_h,n,m)
                di = al_mat[i-1,j-1]+di
                ho = al_mat[i,j-1]
                ve = al_mat[i-1,j]
                al_mat[i,j] = min(di,min(ho,ve))
                p_mat[i,j] = self.Pointer(di,ho,ve)
        if(g.set_VCD == False):
            freq,inten,old_freq = self.Backtrace_IR(p_mat,al_mat,n,m,freq,inten,exp_freq,bond_l,bond_h)
        else:
            freq,inten,inten_VCD,old_freq = self.Backtrace_VCD(p_mat,al_mat,n,m,freq,inten,inten_VCD,exp_freq,bond_l,bond_h)
        returnvalue = al_mat[n-1,m-1]/(n) ##NORMALIZE VALUE BY NUMBER OF THEORETICAL PEAKS. MIGHT BE CHANGED FOR MORE
        if(g.set_VCD==False):
            return returnvalue, old_freq, freq, inten
        else:
            return returnvalue, old_freq, freq, inten, inten_VCD

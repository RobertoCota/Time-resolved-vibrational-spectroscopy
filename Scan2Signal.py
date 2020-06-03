# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 15:52:42 2018

@author: Cota
"""


import numpy as np

from Supporting_Functions import Scans_Average

eps = np.finfo(float).eps




#Subtracts the weighted average spectra of all times before t0 (default -2)
def neg_Substraction(Signal, Error, time, t0):
    
    for i in range (len(time)):
        if (time[i] < t0): t_negat = i 
        
    neg_sig1 = Signal[:,0:t_negat+1,:]
    neg_err1 = Error[:,0:t_negat+1,:]
        
    Weight = 1.0/(neg_err1*neg_err1)
#    Weight = 1.0/(neg_err1*neg_err1)
        
    Weighted_mean = np.average(neg_sig1, axis = 1, weights = Weight) 
        
    Weighted_variance = 1.0/(np.sqrt(np.sum((Weight), axis = 1)))
#    Weighted_variance = np.sqrt(np.average((neg_sig1 - Weighted_mean[:,None,:])**2, weights = Weight,axis=1)) 
            
    for i in range (len(Signal)):
        Signal[i] = Signal[i] - Weighted_mean[i]
        Error[i] = np.sqrt((Error[i]*Error[i]) + (Weighted_variance[i]*Weighted_variance[i]))
        
    return Signal, Error
   
   

  
class Scan2Sig:
    def __init__(self, data, Scans):
        
        if (Scans != 0): 
            N = data.nscans
            A = np.arange(N)
            Delete_scans = np.delete(A, Scans)
            data.s = np.delete(data.s,Delete_scans)
            N = len(data.s)
            print('The following scans were deleted: ', Delete_scans, '\n')
        else:
            N = data.nscans
            print('All the scans are used for the analysis.\n')
            #print (Delete_scans)

        pixel = data.specs.p4v[0]
        wavenumber = data.specs.p4v[1]
        
        self.p = pixel
        self.v = wavenumber

        
        if (N > 0):
            t_check = True
            for i in range (N - 1):
                Comp_t = np.array_equal(data.s[i].t,data.s[i+1].t)  
                if (Comp_t == False): t_check = False
            if (t_check==True): time = data.s[0].t
        else:
            time = data.s[0].t
        
        t_array, E = np.unique(time, return_counts=True)        
        t_pol = np.unique(E)
        if t_pol[0] == 2: time = t_array

        s_par = np.zeros([N,len(time),len(pixel)], dtype = 'float64')
        s_pare = np.zeros([N,len(time),len(pixel)], dtype = 'float64')
        
        s_per = np.zeros([N,len(time),len(pixel)], dtype = 'float64')
        s_pere = np.zeros([N,len(time),len(pixel)], dtype = 'float64')

        for i in range (N):
            for j in range (2*len(time)):
                if(data.s[i].pol[j] == 0): 
                    s_par[i][j] = data.s[i].s2[j]
                    s_pare[i][j] = data.s[i].s2e[j]
                    
                if(data.s[i].pol[j] == 1): 
                    s_per[i][j-len(time)] = data.s[i].s2[j]
                    s_pere[i][j-len(time)] = data.s[i].s2e[j]
        
        
        ###Convert transmittance to transient signal
        S1 = - np.log(s_par)
        S1e = s_pare/s_par

        S2 = - np.log(s_per)
        S2e =  s_pere/s_per
        ###############################
        ###############################
        
        
        ###Delete the extreme delay times to avoid scattering
        t = time[1:len(time)]
        
        S1 = S1[:,1:len(time),:]
        S1e = S1e[:,1:len(time),:]
        S2 = S2[:,1:len(time),:]
        S2e = S2e[:,1:len(time),:]
        ###############################
        ###############################

        ###Substract scatter calculated from measurements 
        ###at t < t0
        S1, S1e = neg_Substraction(S1, S1e, t, t0 = -2.0)
        S2, S2e = neg_Substraction(S2, S2e, t, t0 = -2.0)
        ###############################
        ###############################

        
        # Weighted Scan Average and weighted variance
        S1, S1e = Scans_Average(S1, S1e)
        S2, S2e = Scans_Average(S2, S2e)
        ###############################
        ###############################


        ###Scan for undefined numbers or infinities and covert them
        ###to eps-zero = 2.22e-16
        S1[np.isnan(S1) == True] = eps
        S1[np.isinf(S1) == True] = eps

        S1e[np.isnan(S1e) == True] = eps
        S1e[np.isinf(S1e) == True] = eps

        S2[np.isnan(S2) == True] = eps
        S2[np.isinf(S2) == True] = eps
        
        S2e[np.isnan(S2e) == True] = eps
        S2e[np.isinf(S2e) == True] = eps
        ###############################
        ###############################

        
        self.t = t
        self.S1 = S1
        self.S1e = S1e
        self.S2 = S2
        self.S2e = S2e


###########
#Raw = ReadIris.Read('20180913_Water.dat')
#S = Scan2Sig(Raw, 0)
############

#def Scan_to_Signal(Data, read_scans):
#    Signal = Scan2Sig(Data, read_scans)
#    return Signal





# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 11:20:56 2018

@author: Cota
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.optimize import newton



def v2p(Signal, w):
    idx = (np.abs(Signal.v - w)).argmin()
    return idx
    
def p2v(Signal, pix):
    idx = (np.abs(Signal.p - pix)).argmin()
    return Signal.v[idx]
    
def t2pos(t, y):
    idx = (np.abs(t - y)).argmin()
    return idx    
    
def Swap_pol(S):
    s1 = S.S2
    s1e = S.S2e
    s2 = S.S1
    s2e = S.S1e
    return s1, s1e, s2, s2e
    
def Free_Rot(S):
    Siso = (S.S1 + 2.0*S.S2)/3.0
    Sisoe = np.sqrt((S.S1e*S.S1e) + 4*(S.S2e*S.S2e))/3.0
    return Siso, Sisoe
   
def ShiftT(Signal, root):
    print ('t has been shifted by: %s\n' % (-1.0*root))
    return Signal.t - root
    
def delete_pixel(Signal, p):
    
    Signal.p = np.delete(Signal.p, p, axis = 0)

    Signal.v = np.delete(Signal.v, p, axis = 0)
    
    Signal.S1 = np.delete(Signal.S1, p, axis = 1)
    Signal.S1e = np.delete(Signal.S1e, p, axis = 1)
    
    Signal.S2 = np.delete(Signal.S2, p, axis = 1)
    Signal.S2e = np.delete(Signal.S2e, p, axis = 1)
    
    return Signal.v, Signal.S1, Signal.S1e, Signal.S2, Signal.S2e
   
###Given a matrix of n-scans, m-delay times and l-wavenumbers
###The next function averages over the n-scans and returns
###a matrix with mxl dimensions.
### Weighted Scan Average and weighted variance
def Scans_Average(Signal, Error):
    
    Weight = 1.0/(Error*Error)
    Average_Signal = np.average(Signal, axis = 0, weights = Weight)
    
    Weighted_Variance = 1.0/(np.sqrt(np.sum((Weight), axis = 0)))
#    Weighted_Variance = np.sqrt(np.average((Signal - Average_Signal)**2, weights = Weight,axis=0)) 

    return Average_Signal, Weighted_Variance
    
#def Endt_Spec(Sig, Sig_e, n):
#    Weight = 1.0/(Sig_e[len(Sig)-3:,:]*Sig_e[len(Sig)-3:,:])
#    Signal = np.average(Sig[len(Sig)-3:,:], axis = 0, weights = Weight)
    
    
    
def Plot_Raw(S, W, T, Show=True):
    
    A_iso, A_isoe = Free_Rot(S)
    
    p = v2p(S,W)

    norm = max(np.abs(A_iso[:,p])) 
#    norm = 1.0
    
#    TempP = np.asarray(len(T),dtype = 'int')
#    n=0
#    for i in T:
#            idx = t2pos(S.t,i)
#            TempP[n] = idx
#            n += 1
    
    X = S.t
    YA = S.S1[:,p]/norm
    YB = S.S2[:,p]/norm
    Yiso = A_iso[:,p]/norm
    
    Cs = interp1d(X,Yiso)    
#    print (Cs)
    def f(x):
        return Cs(x) + 0.7
        
    root = newton(f,-0.09)


    if (Show == True):


        #####    
        fig = plt.figure(figsize = (8,5))
        ax1 = fig.add_subplot(111)
        ax1.set_title('Transient absorption at %s cm^1' % S.v[p])
        
        
        ax1.plot(X,YA*norm, '.', color = 'b', label = 'Parallel')
        ax1.plot(X,YB*norm, '.', color = 'g', label = 'Perpendicular')
        ax1.plot(X,Yiso*norm,'.', color = 'r', label = 'Free rotation')
        ax1.plot(X, Cs(X)*norm, color = 'r', label = 'Interpolation') 
        ax1.axvline(x = root, color = 'k')
#        ax1.text(0,0.1,r'$\Delta \tau_0 = %s$' % root, size = 15)
        
        ax1.set_xlabel('Delta Time')
        ax1.set_xlim(-1,50)
#        ax1.set_ylim(Yiso[-1]*norm - 0.0003, Yiso[-1]*norm + 0.0005)
#        ax1.set_xscale('log')
    
        ax1.legend(loc = 4)
        plt.show()


    
        #####    
        fig = plt.figure(figsize = (8,5))
        ax1 = fig.add_subplot(111)
        ax1.set_title('Transient absorption at %s cm^1' % S.v[p])
        
        
        ax1.plot(X,YA*norm, '.', color = 'b', label = 'Parallel')
        ax1.plot(X,YA*norm, color = 'b', label = '')
        ax1.plot(X,YB*norm, '.', color = 'g', label = 'Perpendicular')
        ax1.plot(X,YB*norm, color = 'g', label = '')
        ax1.plot(X,Yiso*norm,'.', color = 'r', label = 'Free rotation')
        ax1.plot(X, Cs(X)*norm, color = 'r', label = 'Interpolation') 
        ax1.axvline(x = root, color = 'k')
#        ax1.text(0,0.1,r'$\Delta \tau_0 = %s$' % root, size = 15)
        
        ax1.set_xlabel('Delta Time')
        ax1.set_xlim(-1,30)
        ax1.set_ylim(Yiso[-1]*norm - 0.004, Yiso[-1]*norm + 0.004)
#        ax1.set_xscale('log')
    
        ax1.legend(loc = 4)
        plt.show()

        
        
        
        ####    
        fig = plt.figure(figsize = (8,5))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
    
        X = S.v
        Y = A_iso[15]
        
        Weight1 = 1.0/(S.S1e[len(S.t)-3:,:]*S.S1e[len(S.t)-3:,:])
        Y1 = np.average(S.S1[len(S.t)-3:,:], axis = 0, weights = Weight1)
        
        Weight2 = 1.0/(S.S2e[len(S.t)-3:,:]*S.S2e[len(S.t)-3:,:])
        Y2 = np.average(S.S2[len(S.t)-3:,:], axis = 0, weights = Weight2)
    
    
        ax1.plot(X,Y1, color = 'b', label = 'Parallel')
        ax1.plot(X,Y2, color = 'g', label = 'Perpendicular')
        
        ax1.set_xlabel('Wavenumbers')
        ax1.set_xlim(3150, 3650)
        ax1.legend(loc = 4)
    
        new_tick_locations = S.v
        
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(S.p)
        ax2.set_xlabel('Pixel')
    
    
        ax2.text(min(X),(max(Y1)+min(Y2))/2.0, 'Average last 3 delay times')
        plt.show()
        
        
        
        
        
        #####    
        fig = plt.figure(figsize = (8,5))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
    
        X = S.v
        
        Y1 = S.S1[2,:]
        Y2 = S.S2[2,:]
    
    
        ax1.plot(X,Y1, color = 'b', label = 'Parallel')
        ax1.plot(X,Y2, color = 'g', label = 'Perpendicular')
        
        ax1.set_xlabel('Wavenumbers')
        ax1.set_xlim(3150, 3650)
        ax1.legend(loc = 4)
    
        new_tick_locations = S.v
        
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(S.p)
        ax2.set_xlabel('Pixel')
    
    
        ax2.text(X[int((len(S.v))/2.0 + 3)],0,'Delay time = - 10 ps')
        plt.show()    
        
        
        
        
        
        
        #####    
        fig = plt.figure(figsize = (8,6))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
    
        X = S.v
    
        colors = plt.cm.jet(np.linspace(0,1,len(T)))
        n = 0
        for i in T:
            idx = t2pos(S.t,i)
#            print (idx)
            Y = A_iso[idx]
            Ye = A_isoe[idx]
            ax1.plot(X, Y, lw = 2, color = colors[len(T) - n - 1], label = '%s' % (round(S.t[idx],2)))
            ax1.errorbar(X, Y, yerr=Ye, color = colors[len(T) - n - 1], fmt='.')
            n += 1
    
        
        ax1.set_xlabel('Wavenumbers')
#        ax1.set_xlim(3100, 3500)
        ax1.set_ylim(-0.09, 0.04)
        ax1.legend(loc = 3, fontsize = 9)
        
        
        new_tick_locations = S.v
        
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(S.p)
        ax2.set_xlabel('Pixel')
        ax1.grid(True)
        plt.show()


    return root
    
    
    
    
def Animate(S):
    
    X = S.v
    
    A_iso, A_isoe = Free_Rot(S)
    
    T = S.t[15:]
    
    
    
    
    n = 0
    for i in T:
        idx = t2pos(S.t,i)
        
        
        fig = plt.figure(figsize = (10,8))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        
        
        colors = plt.cm.jet(np.linspace(0,1,len(T)))
        
        ax1.set_xlabel('Wavenumbers')
        ax1.set_xlim(2390, 2608)
        ax1.set_ylim(-0.06, 0.015)
            
            
        new_tick_locations = S.v
            
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(S.p)
        ax2.set_xlabel('Pixel')
        ax1.grid(True)    
        
        
        Y = A_iso[idx]
        Ye = A_isoe[idx]
        
        ax1.plot(X, Y, color = colors[len(T) - n - 1], lw = 2,label = '%s' % (round(S.t[idx],2)))
        ax1.errorbar(X, Y, yerr=Ye, color = colors[len(T) - n - 1], fmt='.')
        ax1.plot(X, 0.0*X, 'k--', label='')        
        
        
        ax1. text(2410, -0.05, r'$\Delta \tau =$ %s ps' % (round(S.t[idx],2)), fontsize = 18)        
        
        plt.savefig('Frames_NaOH/%s.png' % n)
#        plt.clf()
#        plt.cla()

        n += 1


        
#    ax1.legend(loc = 3, fontsize = 9)
#    plt.show()
    
    

def Decay(t, A, B, C):
    D = A*np.exp(-t/B)+C
    return D

def DecayFix(t, A, B):
    D = A*np.exp(-t/B)
    return D



def Double_Decay(t, A1, tau1, A2, tau2, C):
    D = A1*np.exp(-t/tau1)+ A2*np.exp(-t/tau2) + C
    return D
    
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 12:09:25 2019

@author: Cota
"""

### These modules are needed to read the irisMeet file (Lucy version)
### ReadIris reads and stores all the related information in a variable 
### with several methods and attributes
### Scan2Sig converts the experimentally measured transmitance into
### transient absorption and averages the signal from different scans
### Supporting_Functions contains diverse functions that are needed
### to visualize and simplify the data
from Read_IrisMeet import FileIRISmeet
from Scan2Signal import Scan2Sig
import Supporting_Functions as FUn


### Other modeles needed
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sympy import symbols, Matrix, lambdify, exp
from scipy.optimize import minimize, least_squares, leastsq, curve_fit


#%%
################################
####    Create some lists to import the experimental data
################################ 

cd = 'Data/'

files = ['20190601_NaOD_1000mM.dat',
         '20190601_NaOD_2000mM.dat',
         '20190601_NaOD_3000mM.dat',
         '20190601_NaOD_4500mM.dat',
         '20190601_NaOD_6000mM.dat',
         '20190601_NaOD_8000mM.dat']


####    The following array is created (based on the 
####    experience of the user) to drop the wrong scans
####    Due to temperature variations or air flow in the setup
ScanG = [0,
         0,
         [0,2,3,4,6],
         [0,1,2,3],
         [1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19],
         [0,2,3,4,5,6]]


####    Concentrations 
Sampler = ['1.0 M','2.0 M','3.0 M','4.5 M','6.0 M','8.0 M']


All_Par, All_ParE, All_Per, All_PerE, All_Iso, All_IsoE = [],[],[],[],[],[]
All_Time = []

l=0
for sample in files:
    
    print ('\n')    
    print ('{}\nConcentration: {}\n{}'.format('-'*61,Sampler[l],'-'*20))    
    
    ###Read the file given by the IrisMeet software
    dat = FileIRISmeet(cd + sample)
    
    Sx = Scan2Sig(dat,ScanG[l])       

    ### Swap polarization
    SigA, ErrA, SigB, ErrB = Sx.S1, Sx.S1e, Sx.S2, Sx.S2e 
    Sx.S1, Sx.S1e, Sx.S2, Sx.S2e = SigB, ErrB, SigA, ErrA
    
    ###############################################################################        
    ####################      VISUALIZING THE RAW DATA      #######################
    
    w0 = 3380
    t_plot = np.asarray([-10, 0.45, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 3, 4, 8, 15, 50, 100])
    
    t_nod = FUn.Plot_Raw(Sx, w0, t_plot, Show=False)     #The 4th parameter enables to show the figures: default True

    ###Shift the time to t where A_iso(w0,t)/A_iso^max(w0, t) = - 0.7  
    Sx.t = FUn.ShiftT(Sx,t_nod)
    
    ###Delete pixels 
    Sx.v, Sx.S1, Sx.S1e, Sx.S2, Sx.S2e = FUn.delete_pixel(Sx, [0,1,2,3,4,5,6,7,8,9,25,26,27,28,29,30,31]) 
    
    ###Show raw data
    FUn.Plot_Raw(Sx, w0, t_plot, Show=False)
    
    ###Calculate the isotropic transient absorption
    Ax_iso, Ax_isoe = FUn.Free_Rot(Sx) 
        
    ###############################################################################        
    ########      DEFINING THE THE SPECTRAL AND TIME RANGE ANALYSIS      ########## 
    #########################      FOR THE ANALYSIS      ##########################    
    ###############################################################################        
    
    
    ######      Reduce the spectral range to the window of interest
    f_indx_max = FUn.v2p(Sx,3600)
    f_indx_min = FUn.v2p(Sx,3100)    
    
    freq = Sx.v[f_indx_max:f_indx_min + 1]
    np.savetxt('Matrices/Frequency.dat',freq,fmt='%s', delimiter='\t')    
    
    
    ######      Reduce the time range to the window of interest
    t_indx = FUn.t2pos(Sx.t,0.2)
    t_indx2 = FUn.t2pos(Sx.t,30.5)+1
    
    
    t = Sx.t[t_indx:t_indx2]
       
    np.savetxt('Matrices/TimeRaw_%smM.dat' % int(1000.0*float(Sampler[l].split(' ')[0])), t, fmt='%s', delimiter='\t')
       
    Para, Para_error = Sx.S1[t_indx:t_indx2,f_indx_max:f_indx_min + 1], Sx.S1e[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    Perp, Perp_error = Sx.S2[t_indx:t_indx2,f_indx_max:f_indx_min + 1], Sx.S2e[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    
    # Create the matrix of anisotropy out of parallel and perpendicular
    RRaw = (Para - Perp)/(Para + 2.0*Perp)
    RRaw_Dev = np.sqrt(np.power((3.0*Perp)/np.power(Para+2*Perp,2),2)*(Para_error*Para_error) + np.power((3.0*Para)/np.power(Para+2*Perp,2),2)*(Perp_error*Perp_error))
    
    
    np.savetxt('Matrices/AniRaw_%smM.dat' % int(1000.0*float(Sampler[l].split(' ')[0])), RRaw, fmt='%s', delimiter='\t')
    np.savetxt('Matrices/AniRawStDev_%smM.dat' % int(1000.0*float(Sampler[l].split(' ')[0])), RRaw_Dev, fmt='%s', delimiter='\t')
    
    # In this block, I create and save some matrices to access them in the future
    # to create figures for presentations or publications
    ###############################################################################        
    
    
    
    
    
    ######      Reduce the time range to the window of interest
    t_indx = FUn.t2pos(Sx.t,0.44)
    t_indx2 = FUn.t2pos(Sx.t,30.5)+1
    
    t = Sx.t[t_indx:t_indx2]
    
    np.savetxt('Matrices/Time_%smM.dat' % int(1000.0*float(Sampler[l].split(' ')[0])),t,fmt='%s', delimiter='\t')    
    
    print ('The analysis is set to the spectral range of %s - %s cm^1\nstarting from %s ps\n'% (int(min(freq)),int(max(freq)), round(t[0],2)))
    print('{}'.format('-'*61))
    
    Para, Para_error = Sx.S1[t_indx:t_indx2,f_indx_max:f_indx_min + 1], Sx.S1e[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    Perp, Perp_error = Sx.S2[t_indx:t_indx2,f_indx_max:f_indx_min + 1], Sx.S2e[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    
    Iso_signal = Ax_iso[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    Iso_error = Ax_isoe[t_indx:t_indx2,f_indx_max:f_indx_min + 1]
    
    
    # As in the previous block, here I create and save some matrices to access 
    # them in the future to create figures for presentations or publications
    ###############################################################################    
    
    ###############################################################################   
    # Save all matrices as elements of multiple lists
    All_Par.append(Para)
    All_ParE.append(Para_error)
    
    All_Per.append(Perp)
    All_PerE.append(Perp_error)
    
    All_Iso.append(Iso_signal)
    All_IsoE.append(Iso_error)
    
    All_Time.append(t)
    
    ###############################################################################

    l+=1
    
    ########        END OF READING AND STRUCTURING DATA         ###################
    ###############################################################################
    





#%%
###############################################################################        
###########            SINGULAR VALUE DECOMPOSITION               #############
###############################################################################        

# In some cases a SVD analysis can help to understand the undelying structure
# of the data.

print ('\n')
print ('Singular value decomposition (in %)\n{}'.format('-'*61))

print ('{:11}|{:11}|{:11}|{:11}|{:11}\n{}'.format('Conc', '1st Comp', '2nd Comp', '3rd Comp', '4th Comp', '-'*61))

l=0
for Iso_signal in All_Iso:
    u, s, vh = np.linalg.svd(Iso_signal, full_matrices=True, compute_uv=True)
    print( '{:11}|{:11.2f}|{:11.2f}|{:11.2f}|{:11.2f}'.format(Sampler[l],*((s/sum(s))[:4]*100)))
    l+=1
print ('{}'.format('-'*61))

# >99% of the structure is well represented with the first 3 components
# for all the samples under study here.

###############################################################################        





#%%
###############################################################################        
############        DEFINING A KINETIC MODEL          #######################
###############################################################################     
   
k1, ks = symbols('k1 ks')
TDelay = symbols('t')
A0, B0, C0 = symbols('A0 B0 C0')

KineticModel = Matrix([
                      [ -k1,            0.0,        0.0], 
                      [ +k1,            -ks,        0.0],
                      [ 0.0,            +ks,        0.0]])
                      
# This model represents a unique excited state that decays exponentially to an
# intermidiate state with a decay rate of k1. The system subsequently decays 
# to the (thermal) ground state with a transition rate ks.
# This model descibes well the case of pure water.

P0 = Matrix([A0,B0,C0]) # 

EigenVectors, EigenValues = KineticModel.diagonalize()

Popt = (EigenVectors@exp(EigenValues*TDelay)@(EigenVectors**-1)@P0).T

cascade_model = lambdify((TDelay,k1,ks,A0,B0,C0),Popt,'numpy')



###############################################################################   
##########      A more appropriate model for aqueos ions        ###############
# Not neccesarily the same for all ions but sufficient for the samples under
# study here.

k1, k2, ks = symbols('k1 k2 ks') #Decay rates
TDelay = symbols('t') 
A0, B0, C0, D0 = symbols('A0 B0 C0 D0') #Initial conditions

KineticModel = Matrix([
                      [ -k1,        0.0,        0.0,        0.0], 
                      [ 0.0,        -k2,        0.0,        0.0], 
                      [ +k1,        +k2,        -ks,        0.0],
                      [ 0.0,        0.0,        +ks,        0.0]])

# This model represents two initially excited states that decay exponentially
# with decay rates k1 and k2 to an intermidiate state, which subsequently
# decays to the (thermal) ground state

P0 = Matrix([A0,B0,C0,D0]) # 

EigenVectors, EigenValues = KineticModel.diagonalize()

Popt = (EigenVectors@exp(EigenValues*TDelay)@(EigenVectors**-1)@P0).T

population_model = lambdify((TDelay,k1,k2,ks,A0,B0,C0,D0),Popt,'numpy')


###############################################################################   
###############################################################################   




def ResidualsIso(Parameters, Lambda_Pop, Pop0, freq, t, data, data_e, borders, Fix, Fix_Spec, FullOut):
############################################################################### 
    Rates = Parameters[:3]
###############################################################################    
    Popt = (Lambda_Pop(t,*Rates,*Pop0)[0,:,:]).T
###############################################################################    

    if Fix == True:
        
        Fix_Levels = (Fix_Spec[:,0]).astype('int')
        
        Fix_Spec = (np.delete(Fix_Spec, 0, axis =1)).astype('float64')
        
        if len(Fix_Levels) == len(Pop0):
            Signatures = Fix_Spec  
       
        if len(Fix_Levels) < len(Pop0):

            Specs = np.empty_like(Fix_Spec, dtype = 'float64')        
        
            Specs[0] = Fix_Spec[0]*Parameters[3]
            Specs[1] = Fix_Spec[1]
#            Specs[2] = Fix_Spec[2]

            Signal = Parameters[4:19], Parameters[19:]
#            print (Signal)
                            
            Insert_Spectra = np.delete(np.arange(len(Pop0)), Fix_Levels)
                
            for j in range(len(Insert_Spectra)):
                Specs  = np.insert(Specs, Insert_Spectra[j], Signal[j],axis = 0)
                
            Signatures = Specs 
            
    Constraints = 1000.0*sum(np.absolute(Parameters[:4]-borders[:,0])-(Parameters[:4]-borders[:,0])) + 1000.0*sum(np.absolute(borders[:,1]-Parameters[:4])-(borders[:,1]-Parameters[:4]))

    residual = (np.power((data - Popt@Signatures)/data_e,2.0)).reshape(-1) + Constraints

    if FullOut == False:
        return residual 
        
    if FullOut == True:
        return residual, Popt, Signatures
###############################################################################        
###############################################################################



def ResidualsIsoGlobal(Parameters, GlobalP, Lambda_Pop, Pop0, freq, t, data, data_e, borders, Fix, Fix_Spec, FullOut):
############################################################################### 
    Rates = GlobalP[0], GlobalP[1], Parameters[:1]
###############################################################################    
    Popt = (Lambda_Pop(t,*Rates,*Pop0)[0,:,:]).T
###############################################################################    

    if Fix == True:
        
        Fix_Levels = (Fix_Spec[:,0]).astype('int')
        
        Fix_Spec = (np.delete(Fix_Spec, 0, axis =1)).astype('float64')
        
        if len(Fix_Levels) == len(Pop0):
            Signatures = Fix_Spec  
       
        if len(Fix_Levels) < len(Pop0):

            Specs = np.empty_like(Fix_Spec, dtype = 'float64')        
        
            Specs[0] = Fix_Spec[0]*Parameters[1]
            Specs[1] = Fix_Spec[1]
#            Specs[2] = Fix_Spec[2]

            Signal = Parameters[2:17], Parameters[17:]
#            print (Signal)
                            
            Insert_Spectra = np.delete(np.arange(len(Pop0)), Fix_Levels)
                
            for j in range(len(Insert_Spectra)):
                Specs  = np.insert(Specs, Insert_Spectra[j], Signal[j],axis = 0)
                
            Signatures = Specs 
            
    Constraints = 1000.0*sum(np.absolute(Parameters[:2]-borders[:,0])-(Parameters[:2]-borders[:,0])) + 1000.0*sum(np.absolute(borders[:,1]-Parameters[:2])-(borders[:,1]-Parameters[:2]))

    residual = (np.power((data - Popt@Signatures)/data_e,2.0)).reshape(-1) + Constraints

    if FullOut == False:
        return residual 
        
    if FullOut == True:
        return residual, Popt, Signatures
###############################################################################        
###############################################################################





# Modeling the isotropic transient signal using the two-excited-state model
# the paraments are calculated for each sample separately
#%%
print ('\nModeling the isotropic transient signal:\n{}'.format('-'*65))
print ('{:12}|{:12}|{:12}|{:12}|{:12}\n{}'.format('Conc', '1/k1', '1/k2', '1/ks', 'Residual', '-'*65))

for j in range (len(files)):
    
    t, Iso_signal, Iso_error = All_Time[j], All_Iso[j], All_IsoE[j] 

    # Set initial population array
    N0 = np.array([0.5, 0.5, 0.0, 0.0])
           
    # Initial guess for decay rates and population rate      ############
    Rate0 = np.array([1.0/0.7, 1.0/0.3, 1.0/0.9, 0.5])
    Borders = np.array([[0.0,10.0], [0.0, 10.00], [0.0,10.0], [0.0,10.0]])
    
    
    WSpecs = np.loadtxt('IsoSignals_Water.dat', unpack=True)
    
    
    ###########         Concatenate all parameters and borders       ##############
    Parameters0 = np.concatenate((Rate0, 0.5*WSpecs[0], 0.0*Iso_signal[0]), axis = 0)
    ###############################################################################   
                  
    ##########Fix the spectral shape of the thermally heated ground state##########
    FixSpecs = True
    
    Weights = 1.0/(Iso_error*Iso_error)
    Iso_Heat = np.average(Iso_signal[(len(t) - 4):,:], axis = 0, weights = Weights[(len(t) - 4):,:])
    
    Set_Spectra = np.vstack((Iso_Heat, Iso_Heat))
    
    
    #print (Set_Spectra)
    
    Set_Spectra = np.insert(Set_Spectra, 0, [2,3], axis = 1)
    
    
    Full_OutPut = False
    Vopt = least_squares(ResidualsIso, Parameters0, args=(population_model, N0, freq, t, Iso_signal, Iso_error, Borders, FixSpecs, Set_Spectra, Full_OutPut),  bounds=(-np.inf, np.inf), method='trf', loss='linear', max_nfev = 3000)
    
    Full_OutPut = True
    Res, IsoNt, IsoSigns = ResidualsIso(Vopt['x'], population_model, N0, freq, t, Iso_signal, Iso_error, Borders, FixSpecs, Set_Spectra, Full_OutPut)
    
    #print (1.0/Vopt['x'])
    
    
    Residual = sum(Res)
    Param_err = np.sqrt(Residual*np.diag(np.linalg.inv((Vopt['jac'].T)@Vopt['jac'])))
    
    
    Valores = Vopt['x']
    
#    print ('\nX2 = %s, at c = %s\n' % (round(Residual,2),Sampler,))
    #print (freq.shape)
    
    Rate = Valores[:3]
    
    Scal = Valores[3]
    
    SpecW = Valores[4:]
    SpecWStDev = Param_err[4:]
    
    
    NameRate = ['T1', 'T2', 'T*']
    
    
    PFit = []
    for i in range(len(Rate)):
        PFit.extend([str(round(1.0/Rate[i],2))+'('+ str(int(100.0*round(Param_err[i]/(Rate[i]**2),2)))+')'])
        
    print ('{:12}|{:12}|{:12}|{:12}|{:12}'.format(Sampler[j], PFit[0], PFit[1], PFit[2], round(Residual,2)))

    IsoFit = IsoNt@IsoSigns


    with PdfPages('Plots/IsoFitSpec_%smM.pdf.pdf' % int(1000.0*float(Sampler[j].split(' ')[0]))) as pdf:
        plt.figure(figsize=(10,4))         
        n = [0,1,2,4,6,8,12,15,20,25,32]      
        colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))  
        for i in range(len(n)):
            plt.plot(freq, Iso_signal[n[i],:], '.', color = colors[i])
            plt.plot(freq, IsoFit[n[i],:], lw = 1.5, color = colors[i], label='%s ps' % (round(t[n[i]],2)))
        plt.xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
        plt.ylabel(r'$\Delta\alpha$',size=20)
        #plt.ylim(-0.1,1.1)
        #plt.xlim(0,4)   
        plt.legend(loc=3)
        plt.grid(True)
    
        pdf.savefig(bbox_inches='tight')
        #plt.show()
        plt.close()


    with PdfPages('Plots/IsoFitCurves_%smM.pdf' % int(1000.0*float(Sampler[j].split(' ')[0]))) as pdf:
        plt.figure(figsize=(10,4))
        n = [8,5,3,1]
        colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
        for i in range(len(n)):
            plt.plot(t, -Iso_signal[:,n[i]], '.', color = colors[i], label='%s cm-1' % int(freq[n[i]]))
            plt.plot(t, -IsoFit[:,n[i]], lw = 1.5, color = colors[i])
        plt.text(15,0.015,'[NaOD] = '+Sampler[j], size = 20)
        plt.xlabel(r'$\mathrm{Time}$ $[\mathrm{ps}]$', size = 20)
        plt.ylabel(r'$-\Delta\alpha$',size=20)
        plt.ylim(0.0,0.02)
        plt.xlim(0,32)
        plt.legend(loc=4)
        plt.grid(True)
        pdf.savefig(bbox_inches='tight')
        #plt.show()
        plt.close()



#%%
# Modeling the isotropic transient signal using the two-excited-state model
# the decay rates k1 and k2 are assumed concentration independent and only 
# ks changes with concentration.
print ('\nModeling the isotropic transient signal (Globally):\n{}'.format('-'*65))
print ('{:12}|{:12}|{:12}\n{}'.format('Iter', '1/k1', '1/k2', '-'*45))


def Global_Analysis(PG, All_Time, freq, All_Iso, All_IsoE, Sampler, Show):
    if Show==False:    
        print ('{:12}|{:12.6f}|{:12.6f}'.format(' ', PG[0], PG[1]))
    if Show==True:
        print ('{:12}|{:12.6f}|{:12.6f}\n{}'.format('Final', PG[0], PG[1], '-'*45))
#    print ('Parameters: ', PG)  
    
    ParamT = 1.0/PG
     
    GlobalResidual = []
    
    TstarF = []
    Signatures = []
    
    if Show == True:
        print ('\n{}'.format('-'*45))
        print ('{:12}|{:12}|{:12}'.format('Conc','1/ks','Residual'))
    
    for idx in range(len(All_Time)):
        
        
        Iso_signal = All_Iso[idx]#.astype('float128')
        Iso_error = All_IsoE[idx]#.astype('float128')
        
        t = All_Time[idx]
        
        ###############################################################################        
        ###############################################################################        
        
        
                     
        N0 = np.array([0.5, 0.5, 0.0, 0.0])
        ###############################################################################   
               
        #######     Initial guest for decay rates and population rate      ############
        Rate0 = np.array([1.0/0.9, 0.5])
        Borders = np.array([[0.0,10.0], [0.0,10.0]])
        
        
        WSpecs = np.loadtxt('IsoSignals_Water.dat', unpack=True)
        
        
        ###########         Concatenate all parameters and borders       ##############
        Parameters0 = np.concatenate((Rate0, 0.5*WSpecs[0], 0.0*Iso_signal[0]), axis = 0)
        ###############################################################################   
                      
        ##########Fix the spectral shape of the thermally heated ground state##########
        FixSpecs = True
        
        Weights = 1.0/(Iso_error*Iso_error)
        Iso_Heat = np.average(Iso_signal[(len(t) - 4):,:], axis = 0, weights = Weights[(len(t) - 4):,:])
        
        Set_Spectra = np.vstack((Iso_Heat, Iso_Heat))
        
        
        #print (Set_Spectra)
        
        Set_Spectra = np.insert(Set_Spectra, 0, [2,3], axis = 1)
        
        
        Full_OutPut = False
        Vopt = least_squares(ResidualsIsoGlobal, Parameters0, args=(ParamT, population_model, N0, freq, t, Iso_signal, Iso_error, Borders, FixSpecs, Set_Spectra, Full_OutPut),  bounds=(-np.inf, np.inf), method='trf', loss='linear', max_nfev = 3000)
        
        Full_OutPut = True
        Res, IsoNt, IsoSigns = ResidualsIsoGlobal(Vopt['x'], ParamT, population_model, N0, freq, t, Iso_signal, Iso_error, Borders, FixSpecs, Set_Spectra, Full_OutPut)
        
        #print (1.0/Vopt['x'])
        
        
        Residual = sum(Res)
        Param_err = np.sqrt(Residual*np.diag(np.linalg.inv((Vopt['jac'].T)@Vopt['jac'])))
        
        
        Valores = Vopt['x']
        
        
        Rate = Valores[:1]
        
        
        GlobalResidual.extend(Res)  
        TstarF.extend(1.0/Rate)
        Signatures.append(IsoSigns)
        
    
        
        if Show == True:
            PFit = []
            for i in range(len(Rate)):
                PFit.extend([str(round(1.0/Rate[i],2))+'('+ str(int(100.0*round(Param_err[i]/(Rate[i]**2),2)))+')'])
                
            print ('{:12}|{:12}|{:12}'.format(Sampler[idx], PFit[0], round(Residual,2)))
            
            IsoFit = IsoNt@IsoSigns
    
#            np.savetxt('Matrices/ParExp_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),Para,fmt='%s', delimiter='\t')    
            np.savetxt('Matrices/SpectralDecomposition_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), (np.insert(IsoSigns,0,freq,axis=0)).T, fmt='%s', delimiter='\t')
    
            with PdfPages('Plots/GlobalSpecDecomposition_%smM.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                plt.figure(figsize=(8.5,6))
                LabelLev = [r'$\mathrm{Bulk\/\/\mathbf{OH}}$',r'$\mathrm{\mathbf{OH}} \cdots \mathrm{OD}^-$','Inter',r'$\mathrm{Thermal\/\/bath}$']
                for i in range(len(N0)):
                    if i != 6:
                        plt.plot(freq, 1000*IsoSigns[i,:], lw=3, label = LabelLev[i])
                Tt = np.arange(3230,3475,1)     
                plt.plot(Tt, 0.0*Tt, 'k--', lw=1.5)
                plt.text(3315, -55, r'$\mathrm{[NaOD]} = %s$ $\mathrm{M}$' % (int(1.0*float(Sampler[idx].split(' ')[0]))), size = 35)
                plt.xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 35)
                plt.ylabel(r'$\Delta \alpha \times 10^3$',size=35)
                plt.ylim(-165,20)
                plt.xlim(3240,3465)
                plt.tick_params(labelsize=24)
                plt.legend(loc='upper center',ncol=3,bbox_to_anchor=(0.5, 1.22),fontsize=18.9)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
    
    
            with PdfPages('Plots/GlobalIsoFitSpec_%smM.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                plt.figure(figsize=(10,4))
                n = [0,1,2,4,6,8,12,15,20,25,32]
                colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
                for i in range(len(n)):
                    plt.plot(freq, Iso_signal[n[i],:], '.', color = colors[i])
                    plt.plot(freq, IsoFit[n[i],:], lw = 1.5, color = colors[i], label='%s ps' % (round(t[n[i]],2)))
                plt.xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                plt.ylabel(r'$\Delta\alpha$',size=20)
                plt.legend(loc=3)
                plt.grid(True)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        
        
            with PdfPages('Plots/GlobalIsoFitCurves_%smM.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                plt.figure(figsize=(10,4))
                n = [8,5,3,1]
                colors = plt.cm.gnuplot(np.linspace(0,1,len(n)+1))
                for i in range(len(n)):
                    plt.plot(t, -Iso_signal[:,n[i]], '.', color = colors[i], label='%s cm-1' % int(freq[n[i]]))
                    plt.plot(t, -IsoFit[:,n[i]], lw = 1.5, color = colors[i])
                plt.text(15,0.015,'[NaOD] = '+Sampler[idx], size = 20)
                plt.xlabel(r'$\mathrm{Time}$ $[\mathrm{ps}]$', size = 20)
                plt.ylabel(r'$-\Delta\alpha$',size=20)
                plt.ylim(0.0,0.02)
                plt.xlim(0,32)
                plt.legend(loc=4)
                plt.grid(True)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()


    if Show==False:
        GlobalResidual = np.array(GlobalResidual)
        return GlobalResidual

    if Show==True:
        GlobalResidual = np.array(GlobalResidual)
        Ts = np.array(TstarF)
        Sig = np.array(Signatures)
        return GlobalResidual, Ts, Sig


Parameters0 = [0.76, 0.27]
Vopt = least_squares(Global_Analysis, Parameters0, 
                     args=(All_Time, freq, All_Iso, All_IsoE, Sampler,False),
                     bounds=(-np.inf, np.inf), method='trf', loss='linear', max_nfev = 3000)


Valores = Vopt['x']
Res, Tstar, SigSpec = Global_Analysis(Vopt['x'], All_Time, freq, All_Iso, All_IsoE, Sampler,True)


AllResidual = sum(Res)
GParamE = AllResidual*np.sqrt(np.diag(np.linalg.inv((Vopt['jac'].T)@Vopt['jac'])))

print ('{}\n'.format('-'*65))




#print ('T1 = %s(%s) ps' % (round(Valores[0],2),round(GParamE[0],2)))
#print ('T2 = %s(%s) ps' % (round(Valores[1],2),round(GParamE[1],2)))


###############################################################################
##########          End of the analysis on isotropic signal        ############
###############################################################################







#%%
###############################################################################
##############       START SELF CONSISTENT ANALYSIS         ###################
###############################################################################
# Parallel and perpendicular signals are analyzed together

def ResidualsSelf(Parameters, Lambda_R, Lambda_Pop, TA, Ts, Pop0, freq, t, DataPar, DataPar_e, DataPer, DataPer_e, Spec):
###############################################################################  
    RateA, Rs = 1.0/TA, 1.0/Ts
    
    Popt = (Lambda_Pop(t,*RateA,Rs,*Pop0)[0,:,:]).T

    Rt = Lambda_R(t,*Parameters)[:,0,:].T
    
    PolPar = ((1.0 + 2.0*Rt)*Popt)@Spec
    PolPer = ((1.0 - Rt)*Popt)@Spec
    
    Res = np.vstack((np.power((DataPar - PolPar)/DataPar_e, 2.0), np.power((DataPer - PolPer)/DataPer_e, 2.0))).reshape(-1)# + Constraints
    
    return Res



def Global_SelfAnalysis(Gpar, GlobalT, Tstar, freq, All_Time, All_Par, All_ParE, All_Per, All_PerE, SigSpec, Sampler, PlotA):

#    print ('GParametros: ', Gpar)

    Residuals = []
    
    Conc = [1.0,2.0,3.0,4.5,6.0,8.0]
    
    if PlotA ==True:
        print ('\nModeling the anisotropy:\n{}'.format('-'*65))
#        print('{}'.format('-'*65))
        print('{:12}|{:15}|{:15}|{:15}|{:15}\n{}'.format('Conc', 'tauW', 'tauOD', 'tau0', 'Residual','-'*65))

    for idx in range(len(Sampler)):
    
#        B = symbols('B')
        B, tau1, tau2, tau3 = symbols('B tau1 tau2 tau3')
        TDel = symbols('t')
        
        epsilon = 1e-16
        
        D = Gpar[0]
        
        R0 = B*exp(-TDel/tau1)
        R1 = B*exp(-TDel/tau2)
        R2 = TDel*epsilon + 0.4
        R3 = D*exp(-TDel/tau3)
        
        R = Matrix([R0, R1, R2, R3])
        
        Lam_Ani = lambdify((TDel,B,tau1,tau2,tau3),R,'numpy')
        
        
        #############################################
        
        N0 = np.array([0.5, 0.5, 0.0, 0.0])    
        
        Parameters0 = np.array([0.4,2.6,1.0,4.5])
        
        #############################################
        
        Ts = Tstar[idx]
        
        t = All_Time[idx]
        
        Para = All_Par[idx]
        Para_err = All_ParE[idx]
        
        Perp = All_Per[idx]
        Perp_err = All_PerE[idx]
        
        Sign = SigSpec[idx]
        
        #############################################
        IsoExp = (Para + 2.0*Perp)/3.0
        IsoExp_StDev = np.sqrt((Para_err*Para_err) + 4*(Perp_err*Perp_err))/3.0
        #############################################
        #############################################
        AniExp = (Para - Perp)/(Para + 2.0*Perp)
        AniExp_Dev = np.sqrt(np.power((3.0*Perp)/np.power(Para+2*Perp,2),2)*(Para_err*Para_err) + np.power((3.0*Para)/np.power(Para+2*Perp,2),2)*(Perp_err*Perp_err))
        #############################################
        
        
        
        SelfVopt = least_squares(ResidualsSelf, Parameters0, args=(Lam_Ani, population_model, GlobalT, Ts, N0, freq, t, Para, Para_err, Perp, Perp_err, Sign),  bounds=(0.0, (0.4,np.inf,np.inf,np.inf)), method='trf', loss='linear', max_nfev = 3000)
        
        
        ResSelf = ResidualsSelf(SelfVopt['x'], Lam_Ani, population_model, GlobalT, Ts, N0, freq, t, Para, Para_err, Perp, Perp_err, Sign)
        TResidual = sum(ResSelf)
        
        ParamSelf = SelfVopt['x']
        ParamSelf_err = np.sqrt(np.diag(TResidual*np.linalg.inv((SelfVopt['jac'].T)@SelfVopt['jac'])))
        
        Residuals.extend(ResSelf)
                
#        print ('Hola')
        
        #############################################
        GlobalK = 1.0/GlobalT
        Ks = 1.0/Ts
                
        Ntt = (population_model(t,*GlobalK,Ks,*N0)[0,:,:]).T
        FitIso = Ntt@Sign
        
        R = Lam_Ani(t,*ParamSelf)[:,0,:].T

        FitPar = ((1.0 + 2.0*R)*Ntt)@Sign
        FitPer = ((1.0 - R)*Ntt)@Sign
        
        RFit = (FitPar - FitPer)/(FitPar + 2.0*FitPer)

        #############################################
        ############FIT ANISOTROPY BULK##############

        TFit = np.arange(0.0,8.0,0.05)
        Rani_Bulk = (Lam_Ani(TFit,*ParamSelf)[:,0,:].T)[:,0]
        
        RFitBulk = np.vstack((TFit,Rani_Bulk))
        np.savetxt('Matrices/AniBulkComp_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), RFitBulk.T,fmt='%s', delimiter='\t')    
        


        #############################################

        #############BULK ANISOTROPY#################
#        FitPar_bulk = ((1.0 + 2.0*R[:,0:1])*Ntt[:,0:1])@Sign[0:1,:]
#        FitPer_bulk = ((1.0 - R[:,0:1])*Ntt[:,0:1])@Sign[0:1,:]
        
#        RFit_bulk = (FitPar_bulk - FitPer_bulk)/(FitPar_bulk + 2.0*FitPer_bulk)
        AniBulk = R[:,0]
        #############################################


        ##########COMPLEMENTARY ANISOTROPY###########
        FitPar_Com = ((1.0 + 2.0*R[:,1:])*Ntt[:,1:])@Sign[1:,:]
        FitPer_Com = ((1.0 - R[:,1:])*Ntt[:,1:])@Sign[1:,:]
        
#        RFit_Com = (FitPar_Com - FitPer_Com)/(FitPar_Com + 2.0*FitPer_Com)
#        AniBulk = R[:,0]
        #############################################

        #############BULK EXPERIMENTAL###############
        ParExp_Bulk = Para - FitPar_Com
        PerExp_Bulk = Perp - FitPer_Com
        
        RExp_Bulk = (ParExp_Bulk - PerExp_Bulk)/(ParExp_Bulk + 2.0*PerExp_Bulk)

        #############################################
        #############################################


        
        if PlotA == True:
            print('{:12}|{:15}|{:15}|{:15}|{:15}'.format(Sampler[idx], str(ParamSelf[1].round(2))+'('+str(ParamSelf_err[1].round(2))+')', str(ParamSelf[2].round(2))+'('+str(ParamSelf_err[2].round(2))+')', str(ParamSelf[3].round(2))+'('+str(ParamSelf_err[3].round(2))+')', round(sum(ResSelf),2)))
#            print (ParamSelf,'\t', round(sum(ResSelf),2))
#            print (ParamSelf_err,'\n')
            
            np.savetxt('Matrices/IsoExp_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),IsoExp,fmt='%s', delimiter='\t')    
            np.savetxt('Matrices/IsoExpStDev_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),IsoExp_StDev,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/ParExp_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),Para,fmt='%s', delimiter='\t')    
            np.savetxt('Matrices/ParExpStDev_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),Para_err,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/PerExp_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),Perp,fmt='%s', delimiter='\t')    
            np.savetxt('Matrices/PerExpStDev_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])),Perp_err,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/PopuMat_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), Ntt,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/IsoFit_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), FitIso,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/ParFit_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), FitPar,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/PerFit_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), FitPer,fmt='%s', delimiter='\t')    

            np.savetxt('Matrices/AniFit_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), RFit,fmt='%s', delimiter='\t')    

        
            with PdfPages('Plots/GlobalSelf_DelayCurves_%smM.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                fig, [ax1, ax2, ax3] = plt.subplots(3, 1, sharex='col', sharey='row', figsize=(10,10))
                fig.subplots_adjust(hspace=0.03)
                
                fplot = 8
                ax1.plot(t,-FitPer[:,fplot], color = 'r', lw = 1.5, label='')
                ax1.errorbar(t, -Perp[:,fplot], yerr = Perp_err[:,fplot], color= 'r', fmt='s')
                ax1.plot(t,-FitIso[:,fplot], color = 'y', lw = 2, label='')
                ax1.errorbar(t, -IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='D')
                ax1.plot(t,-FitPar[:,fplot], color = 'b', lw = 1.5, label='')
                ax1.errorbar(t, -Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='o')
                ax1.text(6.0, 0.016, r'$\nu = $ %s cm$^{-1}$' % int(freq[fplot]), fontsize=15)
                ax1.text(11.0, 0.016, r'[NaOD] = %s M' % Conc[idx], fontsize=15)
                ax1.set_ylim(0.0, 0.02)
                ax1.set_xlim(0, 15)
                #ax1.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax1.set_ylabel(r'$-\Delta\alpha$',size=20)
                #ax1.legend(loc=3)
                
                fplot = 5
                ax2.plot(t,-FitPer[:,fplot], color = 'r', lw = 1.5, label='')
                ax2.errorbar(t, -Perp[:,fplot], yerr = Perp_err[:,fplot], color= 'r', fmt='s')
                ax2.plot(t,-FitIso[:,fplot], color = 'y', lw = 2, label='')
                ax2.errorbar(t, -IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='D')
                ax2.plot(t,-FitPar[:,fplot], color = 'b', lw = 1.5, label='')
                ax2.errorbar(t, -Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='o')
                ax2.text(6.0, 0.016, r'$\nu = $ %s cm$^{-1}$' % int(freq[fplot]), fontsize=15)
                ax2.text(11.0, 0.016, r'[NaOD] = %s M' % Conc[idx], fontsize=15)
                ax2.set_ylim(0.0, 0.02)
                ax2.set_xlim(0, 15)
                #ax2.text(3300, 0.02, r'Isotropic', fontsize=15)
                #ax2.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax2.set_ylabel(r'$-\Delta\alpha$',size=20)
                #ax2.legend(loc=3)
                
                fplot = 3
                ax3.plot(t,-FitPer[:,fplot], color = 'r', lw = 1.5, label='')
                ax3.errorbar(t, -Perp[:,fplot], yerr = Perp_err[:,fplot], color= 'r', fmt='s')
                ax3.plot(t,-FitIso[:,fplot], color = 'y', lw = 2, label='')
                ax3.errorbar(t, -IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='D')
                ax3.plot(t,-FitPar[:,fplot], color = 'b', lw = 1.5, label='')
                ax3.errorbar(t, -Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='o')
                ax3.text(6.0, 0.016, r'$\nu = $ %s cm$^{-1}$' % int(freq[fplot]), fontsize=15)
                ax3.text(11.0, 0.016, r'[NaOD] = %s M' % Conc[idx], fontsize=15)
                ax3.set_ylim(0.0, 0.02)
                ax3.set_xlim(0, 15)
                #ax3.text(3300, 0.02, r'Parallel', fontsize=15)
                ax3.set_xlabel(r'$\mathrm{Time}$ $[\mathrm{ps}]$', size = 17)
                ax3.set_ylabel(r'$-\Delta\alpha$',size=20)
                #ax3.legend(loc=3)
                
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        
        
            with PdfPages('Plots/ResidualsGlobalSelf_DelayCurves_%smM.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                fig, [ax1, ax2] = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(10,8))
                fig.subplots_adjust(hspace=0.03)
                fplot = 5
                
                ax1.plot(t,-FitPer[:,fplot], color = 'r', lw = 1.5, label='')
                ax1.errorbar(t, -Perp[:,fplot], yerr = Perp_err[:,fplot], color= 'r', fmt='s')
                ax1.plot(t,-FitIso[:,fplot], color = 'y', lw = 2, label='')
                ax1.errorbar(t, -IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='D')
                ax1.plot(t,-FitPar[:,fplot], color = 'b', lw = 1.5, label='')
                ax1.errorbar(t, -Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='o')
                ax1.text(6.0, 0.016, r'$\nu = $ %s cm$^{-1}$' % int(freq[fplot]), fontsize=15)
                ax1.text(11.0, 0.016, r'[NaOD] = %s M' % Conc[idx], fontsize=15)
                ax1.set_ylim(0.0, 0.02)
                ax1.set_xlim(0, 15)
                #ax1.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax1.set_ylabel(r'$-\Delta\alpha$',size=20)
                #ax1.legend(loc=3)

                ax2.errorbar(t, FitPar[:,fplot] - Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='.')
                ax2.errorbar(t, FitIso[:,fplot] - IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='.')
                ax2.errorbar(t, FitPer[:,fplot] - Perp[:,fplot] , yerr = Perp_err[:,fplot], color= 'r', fmt='.')
                ax2.plot(t,0.0*t, 'k--', lw = 1, label='')
                #ax2.set_ylim(0.0, 0.02)
                ax2.set_xlim(0, 15)
                #ax2.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax2.set_ylabel(r'$-\Delta\alpha$',size=20)
                #ax2.legend(loc=3)
                
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        
        
            with PdfPages('Plots/GlobalSelf_DelayCurves_%smM_Single.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                plt.style.use('classic')
                fig, ax1 = plt.subplots(1, 1, figsize=(7,5))
                #fig.subplots_adjust(hspace=0.03)
                
                fplot = 5
                    
                ax1.plot(t,-FitPar[:,fplot], color = 'b', lw = 1.5, label='Parallel')
                ax1.errorbar(t, -Para[:,fplot], yerr = Para_err[:,fplot], color= 'b', fmt='o')
                ax1.plot(t,-FitIso[:,fplot], color = 'y', lw = 2, label='Isotropic')
                ax1.errorbar(t, -IsoExp[:,fplot], yerr = IsoExp_StDev[:,fplot], color= 'y', fmt='D')
                ax1.plot(t,-FitPer[:,fplot], color = 'r', lw = 1.5, label='Perpendicular')
                ax1.errorbar(t, -Perp[:,fplot], yerr = Perp_err[:,fplot], color= 'r', fmt='s')
                ax1.text(7.0, 0.010, r'$\nu = $ %s cm$^{-1}$' % int(freq[fplot]), fontsize=15)
                ax1.text(7.0, 0.012, r'[NaOD] = %s M' % Conc[idx], fontsize=15)
                ax1.set_ylim(0.0, 0.02)
                ax1.set_xlim(0, 15)
                #ax1.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax1.set_ylabel(r'$-\Delta\alpha$',size=20)
                ax1.set_xlabel(r'$\mathrm{Delay}$ $\mathrm{time}$ $\mathrm{[ps]}$',size=20)
                ax1.legend(loc=1)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        
        
            with PdfPages('Plots/GlobalSelf_AnisotropyCurves_%smM_Single.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                plt.style.use('classic')
                fig, ax1 = plt.subplots(1, 1, figsize=(7,5))
                #fig.subplots_adjust(hspace=0.03)
                
                fplot = 5
                    
                #AniFit = (FitPar -  FitPer)/(FitPar +  2.0*FitPer)  
                #AniExp = (Para - Perp)/(Para + 2.0*Perp)
                #AniExp_Dev = np.sqrt(np.power((3.0*Perp)/np.power(Para+2*Perp,2),2)*(Para_err*Para_err) + np.power((3.0*Para)/np.power(Para+2*Perp,2),2)*(Perp_err*Perp_err))

                ax1.plot(t,RFit[:,fplot], color = 'b', lw = 1.5, label='Anisotropy')
                ax1.errorbar(t, AniExp[:,fplot], yerr = AniExp_Dev[:,fplot], color= 'b', fmt='o')
                ax1.text(1, 0.0, '%s mM' % int(1000.0*float(Sampler[idx].split(' ')[0])))
                ax1.set_ylim(-0.02, 0.4)
                ax1.set_xlim(0, 7)
                #ax1.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax1.set_ylabel(r'$R(t)$',size=20)
                ax1.set_xlabel(r'$\mathrm{Delay}$ $\mathrm{time}$ $\mathrm{[ps]}$',size=20)
                ax1.legend(loc=1)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        


            with PdfPages('Plots/GlobalSelf_AnisotropyCurves_%smM_Check.pdf' % int(1000.0*float(Sampler[idx].split(' ')[0]))) as pdf:
                #plt.style.use('classic')
                fig, ax1 = plt.subplots(1, 1, figsize=(7,5))
                #fig.subplots_adjust(hspace=0.03)
                
                fplot = 4
                    
                #AniFit = (FitPar -  FitPer)/(FitPar +  2.0*FitPer)  
                #AniExp = (Para - Perp)/(Para + 2.0*Perp)
                AniExp_DevBulk = np.sqrt(np.power((3.0*PerExp_Bulk)/np.power(ParExp_Bulk+2.0*PerExp_Bulk,2),2)*(Para_err*Para_err) + np.power((3.0*ParExp_Bulk)/np.power(ParExp_Bulk+2.0*PerExp_Bulk,2),2)*(Perp_err*Perp_err))                
                ax1.plot(t,AniBulk, color = 'r', lw = 2.8, label='Anisotropy')
                #ax1.plot(t,RFit_bulk[:,fplot], color = 'b', lw = 1.4, label='Anisotropy')
                ax1.errorbar(t, RExp_Bulk[:,fplot], yerr = AniExp_DevBulk[:,fplot], color= 'b', fmt='o')
                #ax1.plot(t,RExp_Bulk[:,fplot], 's', color = 'm', lw = 1.4, label='Anisotropy')
                Ani_Bulk = np.vstack((t,RExp_Bulk[:,fplot],AniExp_DevBulk[:,fplot],AniBulk))                
                
                np.savetxt('Matrices/AniBulk_%smM.dat' % int(1000.0*float(Sampler[idx].split(' ')[0])), Ani_Bulk.T,fmt='%s', delimiter='\t')    
                    
                ax1.text(1, 0.0, '%s mM' % int(1000.0*float(Sampler[idx].split(' ')[0])))
                ax1.set_ylim(-0.02, 0.4)
                ax1.set_xlim(0, 3.5)
                #ax1.set_xlabel(r'$\mathrm{Wavenumbers}$ $[\mathrm{cm}^{-1}]$', size = 20)
                ax1.set_ylabel(r'$R(t)$',size=20)
                ax1.set_xlabel(r'$\mathrm{Delay}$ $\mathrm{time}$ $\mathrm{[ps]}$',size=20)
                ax1.legend(loc=1)
                pdf.savefig(bbox_inches='tight')
                #plt.show()
                plt.close()
        
        
        
    return sum(Residuals)


Param0 = [0.11]

Plot = False
GlobalSelfOpt = least_squares(Global_SelfAnalysis, Param0, args=(Valores, Tstar, freq, All_Time, All_Par, All_ParE, All_Per, All_PerE, SigSpec, Sampler, Plot),  bounds=(0.0, 0.4), method='trf', loss='linear', max_nfev = 3000)

Param = GlobalSelfOpt['x']

Plot = True
Res = Global_SelfAnalysis(Param, Valores, Tstar, freq, All_Time, All_Par, All_ParE, All_Per, All_PerE, SigSpec, Sampler, Plot)


# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:30:22 2018

@author: Cota
"""

import numpy as np

from datetime import date, time


###This function looks for a string "String" given by the user
###in the file "dat" line by line and return the number of line
###where the string was found. Note: Once thestring is found the 
###function stops. So, it is not possible to find in the string 
###is repeated in the file.
def FindinLine(String, dat):
    l = 0 
    Lines = open(dat, 'r').read()
    Lines = Lines.split('\n')
    for line in Lines:
        #print(line)
        if String in line:
            break
        l += 1
    return l
     
###It creates the subclasses for the zero correction parameters
###It returns three classes for the signal1, signal2 and reference data
class ZEROCorr:
    def __init__(self, Sig1, Sig2, Ref):
        self.s1 = Sig1
        self.s2 = Sig2
        self.ref = Ref

class Specifications:
    def __init__(self, Spectro, Grating, CentralP, CentralV, P4V):
        self.SType = Spectro
        self.GLines = Grating
        self.centralp = CentralP
        self.centralv = CentralV
        self.p4v = P4V
        
        #p4v = np.zeros((2,n), )
    
class Spectrum:
    def __init__(self, T, S1, S1e, R1, A1, POL, S2, S2e, R2, A2, S12, S12e, R12, A12):
        self.t = T
        self.s1 = S1
        self.s1e = S1e
        self.ref1 = R1
        self.a1 = A1

        self.pol = POL

        self.s2 = S2
        self.s2e = S2e
        self.ref2 = R2
        self.a2 = A2

        self.s12 = S12
        self.s12e = S12e
        self.ref12 = R12
        self.a12 = A12





                     
class FileIRISmeet:
    def __init__(self, FileIRIS):
        
        ###Opens a file given by the user and store it in the FileComplete
        ###variable to be read and analyze later.
        FileComplete = open(FileIRIS, 'r').read()
        
        
        ###It splits the file in lines and store each line as an string in 
        ###a lines-size list to be read later.
        File = FileComplete.split('\n')
        
        
        ###Read file lines looking for the original file name. It returns
        ###a string, and store it in the 'filename' class
        line = FindinLine('Original filename', FileIRIS)
        LineaFile = File[line]
        self.filename = LineaFile[21:]
        
        
        ###Read file lines looking for general information. It returns
        ###three classes:
        ### -self.data: type 'datetime.date'
        ### -self.time: type 'datetime.time'
        ### -self.version: type 'string'
        line = FindinLine('Delay scans recorded by', FileIRIS)
        Linea1 = File[line]
        Linea1 = Linea1.split()[1:]
        
        dateArray = Linea1[10].split('-')
        dateArray = np.asarray(list(map(int,dateArray)))
        self.date = date(dateArray[2],dateArray[1],dateArray[0])

        timeArray = Linea1[11][:8]
        timeArray = timeArray.split(':')
        timeArray = np.asarray(list(map(int,timeArray)))
        self.time = time(timeArray[0], timeArray[1], timeArray[2])

        self.version = Linea1[4] + ' ' + Linea1[5]
        

        ###Read comments given by the user in the measurements. It returns 
        ###a string with the comments that the user may have insert prior 
        ###measurements. It creates a callable class: 'self.comments'
        line = FindinLine('Comment', FileIRIS)
        Comments = File[line]
        self.comments = Comments[2:]


        ###Number of detector pixels. It returns an integer with the number
        ###of available pixels in the detector. It creates a callable class:
        ###'self.NPixels'
        line = FindinLine('Number of detector pixels used', FileIRIS)
        numberPixels = File[line]
        numberPixels = int(numberPixels[34:])
        self.NPixels = numberPixels
        
        
        ###Pixels recorded in the detector. It returns an array with 1's
        ###and 0's. Each element of the array refers to each pixel. 1 means 
        ###that the i-pixel was active during measurements. 0 meas that the 
        ###i-pixel was turn down during measurements.
        line = FindinLine('SHOWPIX_FLAGS', FileIRIS)
        PixelsInUse = File[line]
        PixelsInUse = PixelsInUse[16:].split(',')
        PixelsInUse = np.asarray(list(map(int,PixelsInUse)))
        self.UsePixels = PixelsInUse


        ###Zero corrections Sig1, Sig2, Ref. It returns three arrays with
        ###the different signas, respectively. Note: the arrays are store 
        ###as subclasses of the 'self.zerocorr' class. 
        line = FindinLine('Zero corrections: ', FileIRIS)

        sig1 = np.zeros(numberPixels, dtype = 'float64')
        sig2 = np.zeros(numberPixels, dtype = 'float64')
        ref = np.zeros(numberPixels, dtype = 'float64')
        for i in range(numberPixels):
            ZeroCorr = File[line + 1 + i][2:].replace(';', '').split()
            sig1[i], sig2[i], ref[i] = float(ZeroCorr[3]), float(ZeroCorr[5]), float(ZeroCorr[7])
        self.zerocorr = ZEROCorr(sig1,sig2,ref)    
        
        
        ###Scan time position and shots per position. It returns two arrays
        ###which contain all the time delay positions and the amount of shots 
        ###that were made at each time delay position, respectively. The 
        ###arrays are stored as classes of the structure: 'self.Delta_t' 
        ###and 'self.NShots'
        Text = FileComplete.split('# # Columns:')
        Text = Text[0].split('Scan positions/#shots: ')
        Text = Text[1].replace('\n','')
        Text = Text.replace('#','')
        Text = Text.replace(',','').split()
        TimePosition = np.zeros(len(Text), dtype = 'float64') 
        ShotsPerPosition = np.zeros(len(Text), dtype = 'int') 
        for i in range (len(Text)):
             TimePosition[i], ShotsPerPosition[i] = Text[i].split(':')[0], Text[i].split(':')[1]
        
        self.Delta_t = TimePosition
        self.NShots = ShotsPerPosition
        
        
        ###Smoothing method: It returns a string, stored at the 'self.smooth'
        ###class
        line = FindinLine('# Interference smoothing: ', FileIRIS)
        Smoothing = File[line][26:]
        self.smooth = Smoothing
        
        
        ###This function creates a class with several subclasses about the 
        ###experiment and spectrograph specifications: 
        ### -self.specs: it cointains 5 different specifications
        ###     -self.SType: Type of spectrograph
        ###     -self.GLines: Lines per mm
        ###     -self.centralp: chosen pixel in the IRIS software
        ###     -self.centralv: wavenumber where the spectrograph is centre
        ###     -self.p4v: It contains an 2xNPixels array. With a relation 
        ###                between pixels and wavenumbers. It is important 
        ###                for further analises.         
        line = FindinLine('# Type of spectrograph:', FileIRIS)
        Spectrograph = File[line].replace('# Type of spectrograph:','').split()[0]
        
        line = FindinLine('# Grating used has', FileIRIS)
        Grating = int(File[line].replace('# Grating used has','').split()[0])
                
        line = FindinLine('# Pixel number in centre of image plane:', FileIRIS)
        CentralPixel = int(File[line].replace('# Pixel number in centre of image plane:','').split()[0])
    
        line = FindinLine('# Spectrograph set at', FileIRIS)
        CentralWavenumber = float(File[line].replace('# Spectrograph set at','').split()[0])
    
        line = FindinLine('# Wavelengths for pixels', FileIRIS)
        pixel4wavenumber = File[line].replace('# Wavelengths for pixels','').split()
        
        pixel = np.zeros(numberPixels,dtype = 'int')
        wavenumber = np.zeros(numberPixels,dtype = 'float64')
        for i in range(numberPixels):
            pixel[i], wavenumber[i] = int(pixel4wavenumber[i].split(':')[0]), float(pixel4wavenumber[i].split(':')[1])        
        p4v = np.ndarray((2,numberPixels), dtype = object)
        p4v[0], p4v[1] = pixel[::1], wavenumber[::1]
        
        self.specs = Specifications(Spectrograph, Grating,  CentralPixel, CentralWavenumber, p4v)




        ###The next function separate the file in the different scans and 
        ###delete the info which has been extracted above. This function 
        ###also test whether the last recorded scan is, in comparison with 
        ###previos scans, completed. If not, the last scan is deleted. All 
        ###the fully-taken scans are stored in 'data_file'
        data_file = FileComplete.split('# scan')
        nscans = len(data_file[1:])
        
        if (len(data_file[nscans].split('\n')) != len(data_file[nscans - 1].split('\n'))):
            data_file = data_file[1:nscans]
            nscans = nscans - 1

        else:
            data_file = data_file[1:nscans+1]
        
        self.nscans = nscans        
        
        print ('\nNumber of complete scans: %d\n' % nscans)
        
        scans = np.ndarray(nscans, dtype = object)
        
        for l in range(nscans):
        
            Datos_file = data_file[l].split('\n')#[1:469]
            
            Datos_file = Datos_file[1:int(len(Datos_file)-2)]
            
            t = np.zeros(int(len(Datos_file)/3), dtype = 'float64')    
            pol = np.zeros(int(len(Datos_file)/3), dtype = 'int')    
    
            s1, s1e, r1, a1 = [], [], [], []
            s2, s2e, r2, a2 = [], [], [], []
            s12, s12e, r12, a12 = [], [], [], []
    
    
            for i in range (len(Datos_file)):
                
                S = list(map(float,Datos_file[i].split()[2:130]))
               
                if 'c0' in Datos_file[i].split(): 
                    t[int(i/3)] = Datos_file[i].split()[1]
                    pol[int(i/3)] = Datos_file[i].split()[130]
    
                    s1.append(S[0::4][::1])
                    s1e.append(S[1::4][::1])
                    r1.append(S[2::4][::1])
                    a1.append(S[3::4][::1])
    
                if 'c1' in Datos_file[i].split(): 
                    s2.append(S[0::4][::1])
                    s2e.append(S[1::4][::1])
                    r2.append(S[2::4][::1])
                    a2.append(S[3::4][::1])
    
                if 'c2' in Datos_file[i].split(): 
                    s12.append(S[0::4][::1])
                    s12e.append(S[1::4][::1])
                    r12.append(S[2::4][::1])
                    a12.append(S[3::4][::1])
    
    
            s1 = np.asarray(s1, dtype = 'float64')  
            s1e = np.asarray(s1e, dtype = 'float64')        
            r1 = np.asarray(r1, dtype = 'float64')  
            a1 = np.asarray(a1, dtype = 'float64')  
            
            s2 = np.asarray(s2, dtype = 'float64')  
            s2e = np.asarray(s2e, dtype = 'float64')        
            r2 = np.asarray(r2, dtype = 'float64')  
            a2 = np.asarray(a2, dtype = 'float64')  
    
            s12 = np.asarray(s12, dtype = 'float64')  
            s12e = np.asarray(s12e, dtype = 'float64')        
            r12 = np.asarray(r12, dtype = 'float64')  
            a12 = np.asarray(a12, dtype = 'float64')  
            
            scans[l] = Spectrum(t, s1, s1e, r1, a1, pol, s2, s2e, r2, a2, s12, s12e, r12, a12)
        
        self.s = scans


#def Read(file):      
#    Data = FileIRISmeet(file)
#    return Data


########
## Create a movie to observe how the transient
## signal changes with the time
########
#with PdfPages('T1.pdf') as pdf:
#    plt.figure()
#    plt.xlim([2350, 2570])
#    plt.ylabel(r'$Absorbance$')
#    plt.xlabel(r'$Wavenumbers [cm^{-1}]$')
#    plt.ylim([0.96, 1.10])
#         
#	#plt.text(0.55,30.0,r'$\kappa_{[HCl]} / \kappa_{[DCl]} \approx \/ \sqrt{2} $')
#
#	#plt.gca().invert_xaxis()
#    for i in range(int(len(Data.s[3].t))):
#        plt.plot(Data.specs.p4v[1], -np.log(Data.s[2].s2[i]),color = 'g',label='')
#        plt.ylim([-0.02, 0.005])
#        plt.savefig('Frames/%s.png' % i)
#        plt.clf()
#        plt.cla()
#
#
#
#	#plt.legend(loc=2)
#    pdf.savefig()
#    plt.close()

#ffmpeg -r 10 -i %d.png -pix_fmt yuv420p outputC.mp4


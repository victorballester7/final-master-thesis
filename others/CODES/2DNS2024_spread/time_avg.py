import scipy
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, sys
from pylab import *
from scipy.interpolate import interp1d


def readtxt(inputfilename,nx,ny):
    f = open(inputfilename,'r')
    field = np.fromfile(f,dtype = 'd',count=int(nx*ny))
    f.close()
    return field

def readslice(inputfilename,nx,ny):
    f = open(inputfilename,'rb')
    f.seek(4)
    field = np.fromfile(f,dtype = 'd',count=int(nx*ny))
    f.close()
    return field


filename1 = 'kspectrum.' 
filename2 = 'Structure_k.'  
filename3 = 'Structure_p.'  
filename4 = 'vectrans_euu.'  
filename5 = 'vectrans_vrt.'

outnum_sta=200
outnum_end=1000
NL=int(512/2)

EN=np.zeros((2,NL))
SK=np.zeros((9,NL))
SP=np.zeros((9,NL))
EF=np.zeros((2,NL))
VF=np.zeros((2,NL))
for line in range(NL):
    EN[0,line] = line+1
    SK[0,line] = (line+1)*2*3.14159/(2*(NL-1))
    SP[0,line] = (line+1)*2*3.14159/(2*(NL-1))
    EF[0,line] = line+1
    VF[0,line] = line+1

for ifile in range(outnum_sta,outnum_end):
    data1 = np.loadtxt(filename1+str("%04d" % ifile)+'.txt')
    data2 = np.loadtxt(filename2+str("%04d" % ifile)+'.txt')
    data3 = np.loadtxt(filename3+str("%04d" % ifile)+'.txt')
    data4 = np.loadtxt(filename4+str("%04d" % ifile)+'.txt')
    data5 = np.loadtxt(filename5+str("%04d" % ifile)+'.txt')
    for line in range(NL):
        EN[1,line]  = EN[1,line]+data1[line]
        for col in range(8):
            SK[col+1,line]  = SK[col+1,line]+data2[line,col]
            SP[col+1,line]  = SP[col+1,line]+data3[line,col]
        EF[1,line]  = EF[1,line]+data4[line,1]
        VF[1,line]  = VF[1,line]+data5[line,1]

for line in range(NL):
    print(line)
    EN[1,line] = (EN[1,line])/(outnum_end-outnum_sta+1)
    EF[1,line] = (EF[1,line])/(outnum_end-outnum_sta+1)
    VF[1,line] = (VF[1,line])/(outnum_end-outnum_sta+1)
    for col in range(8):
            SK[col+1,line]=SK[col+1,line]/(outnum_end-outnum_sta+1)
            SP[col+1,line]=SP[col+1,line]/(outnum_end-outnum_sta+1)


######################################
figure(1)
plt.ylabel("$S_2,S_4,S_6$")
plt.xlabel("r")
plt.xscale("log")
plt.yscale("log")
plt.ylim([0.0001,10000])
plt.plot(SP[0,:],SP[2,:],'b')
plt.plot(SP[0,:],SP[4,:],'r')
plt.plot(SP[0,:],SP[6,:],'k')
#plt.plot(SP[0,:],(SP[0,:]/SP[0,0])**0.666/100,'--k')
#plt.plot(SP[0,:],(SP[0,:]/SP[0,0])**1.333/100,'--k')
#plt.plot(SP[0,:],(SP[0,:]/SP[0,0])**2.000/100,'--k')

######################################
figure(2)
plt.xscale("log")
plt.yscale("log")
plt.ylabel("$S_3,S_5,S_7$")
plt.xlabel("r")
plt.ylim([0.01,10000])
plt.plot(SP[0,:],SP[3,:],'b')
plt.plot(SP[0,:],SP[5,:],'r')
plt.plot(SP[0,:],SP[7,:],'k')
#plt.plot(SP[0,:],(SP[0,:]/SP[0,0])**1.000/100,'--k')
#plt.plot(SP[0,:],(SP[0,:]/SP[0,0])**1.666/100,'--k')

######################################
figure(3)
plt.xscale("log")
plt.yscale("log")
plt.ylim([0.01,100])
plt.ylabel("$S_2,S_4^{1/2},S_6^{1/3}$")
plt.xlabel("r")
plt.plot(SP[0,:],SP[2,:],'b')
plt.plot(SP[0,:],SP[4,:]**0.5,'r')
plt.plot(SP[0,:],SP[6,:]**0.333,'k')
plt.plot(SP[0,:],200*(SP[0,:]/SP[0,0])**0.666/100,'--k')

######################################
figure(4)
plt.xscale("log")
plt.yscale("log")
plt.ylim([0.01,100])
plt.ylabel("$S_4^{1/2}/S_2,S_6^{1/3}/S_2$")
plt.xlabel("r")
plt.plot(SP[0,:],SP[4,:]**0.5/SP[2,:],'r')
plt.plot(SP[0,:],SP[6,:]**0.3333/SP[2,:],'k')

######################################
figure(5)
plt.xscale("log")
plt.yscale("log")
plt.ylim([0.01,100])
plt.ylabel("$S_5,S_7$")
plt.xlabel("S_3")
plt.plot(SP[3,:],SP[5,:]**0.5/SP[2,:],'r')
plt.plot(SP[3,:],SP[7,:]**0.3333/SP[2,:],'k')



plt.show()

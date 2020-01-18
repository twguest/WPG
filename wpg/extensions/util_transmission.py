#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:16:21 2019

@author: -
"""

import imageio
import xrt.backends.raycing.materials as rm
import numpy as np
import numpy as np
from numpy import argmin
import matplotlib.pyplot as plt
import imageio

class kinoform():
    
    def __init__(self,
                 diameter=100.e-6, 
                 Nzones = 50, 
                 delta = 0.00114228425,
                 beta = 0.00091671315, 
                 E = 1000, #eV
                 Npoints=5000):
        
        self.diameter = diameter
        self.Nzones = Nzones
        self.delta = delta
        self.beta = beta
        self.E = E
        self.wl = self.wl(E)
        self.Npoints = Npoints
        self.f = self.diameter*(self.diameter/(4*self.Nzones)) / self.wl
        self.x = np.linspace(-diameter/2.0, diameter/2.0, num=Npoints)
        
        self.m = 2 # an integer multiple m of 2pi phase shift occus at each zone
        
    def wl(self,E):
        #return wavelength for energy E in keV
        return  12.39e-10/E
        
    def xj0(self,j):
        # returns start position of jth zone.  This is essentially the Zone Plate law    
        xval = np.sqrt( 2*j*(self.m*self.wl)*self.f + j**2*(self.m*self.wl)**2 )
        return xval
    
    def xj(self,j):
        # return list of x points spanning the jth zone
        return np.linspace(self.xj0(j),self.xj0(j+1), self.xwpts(j), endpoint=True)
    
    def xw(self,j):
        # return width of jth zone
        return self.xj0(j+1) - self.xj0(j)
    
    def xwpts(self, j):
        # return number of points across jth zone
        return int(self.Npoints * self.xw(j)/self.diameter)

    def zj(self,j):
        # thickness profile of the jth zone
        #  Note that this is computed as the OPD from edge of zone.  
        z = (np.sqrt(self.xj(j)**2+self.f**2)-np.sqrt(self.xj0(j)**2+self.f**2) )/ self.delta
        z = z / self.m  #added for consistency with xj - please check!
        return z
    
    def z(self):
        x, zz = [],  []
        for j in range(self.Nzones):
            x.extend(self.xj(j))
            zz.extend( self.zj(j)  )    
        return zz, x    
    
    def z_symmetric(self):        
        z, x = self.z() 
        return np.concatenate((z[::-1], z)), np.concatenate((x[::-1],x))
               
    def zEquiv(self):
        # return the kinoform profile, without each zone shifted to minimise aborption      
        return (np.sqrt(self.x**2+self.f**2) -self.f) / self.delta  
    
    def phase(self):
        #return phase shift through thickness z, relative to vacuum
        x, p = [],  []
        for j in range (self.Nzones):
            x.extend(self.xj(j))
            p.extend( (2*np.pi*self.delta/self.wl)*self.zj(j))
        
        return p, x 
        
    
def test(diameter=15.0e-6, NumberOfZones=36,pixelSize=50e-9):
    resolution = pixelSize
    NX = int(diameter/resolution)
    kino = kinoform(diameter=diameter, 
                     Nzones = NumberOfZones, 
                     delta = 0.00314069726,  # for Ni at 600 eV (0.0031406)
                     beta = 0.000591546006,   # for Ni at 1000 eV
                     E = 0.600, #eV
                     Npoints=NX)

    fig, ax = plt.subplots()
    scale = 1
    ax.plot(kino.x, kino.zEquiv()/scale, 
            label='Equivalent profile / {}'.format(kino.Npoints))
    zz, xx = kino.z()
    ax.plot(xx, zz, label='Kinoform')
    ax.set(xlabel='x [m]', ylabel='z [m]')
    plt.legend()
    plt.show()
    
    


def compare():
    diameter=15.0e-6
    kino = kinoform(diameter=diameter, 
                     Nzones = 36, 
                     delta = 0.00314069726,  # for Ni at 600 eV (0.0031406)
                     beta = 0.000591546006,   # for Ni at 1000 eV
                     E = 0.600, #eV
                     Npoints=200000)
    
    P = np.loadtxt(open('AllZones.csv', "rb"), delimiter=",", skiprows=1)
    xc, zc = P[:,0]*10e-9*1e6+2.61, P[:,1]
    phase, xp = kino.phase()
    zz, xx = kino.z()
    
    fig, ax = plt.subplots()
    #ax.plot(kino.x, np.multiply(kino.zEquiv(),1e6), label='Equivalent profile / {}'.format(kino.Npoints))
    ax.plot(np.multiply(xx,1.e6), np.multiply(zz,1e6), label='Kinoform', linewidth=0.1)
    ax.plot(xc, np.multiply(zc,1.e6), label = 'other source', linewidth=0.1)
    #ax2 = ax.twinx()
    #ax2.plot(np.multiply(xp,1.e6),phase, label='Kinoform phase shift', linewidth=0.1)
    ax.set(xlabel='x [um]', ylabel='z [um]')

    plt.legend()
    plt.show()
    
    
    fig, ax = plt.subplots()
    ax.plot(np.multiply(xp,1.e6),phase, label='Kinoform phase shift', linewidth=0.1)
    ax.set(xlabel='x [um]', ylabel='phase shift [rad.]')
    plt.legend()
    plt.show()
    

def make2D(profile,dimensions):   
   return np.tile(profile, dimensions)

def nearestValueIndex(x,L):
    # return the index of the value in L that is closest to x
    idx = (np.abs(np.subtract(L,x))).argmin()
    return idx #, L[idx]   

def radialSymmetric(profile,x,crop=False):
    #  
    rdist = distanceFromMedian2D(x)
    z = [profile[nearestValueIndex(val,x)] for val in np.nditer(rdist)]
    z = np.reshape(z,[len(x),len(x)])
    
    if crop is True:
        z[np.where(rdist>np.max(x)/2)] = 0
    
    return z, rdist
    
    
def kinoform2D(diameter=15.0e-6, 
               NumberOfZones=36,
               delta = 0.00189737359,  # for Si3N4 at 600 eV (0.330977 um attenuation length)
               beta =  0.000496833818,   #  for Si3N4 at 600 eV 
               E = 0.600, #eV
               pixelSize=50e-9,
               savePath=None, display=False):

    resolution = pixelSize
    NX = int(diameter/resolution)*2
    print('Pixel size: {}'.format(resolution))
    kino = kinoform(diameter=diameter, 
                     Nzones = NumberOfZones, 
                     delta = delta,
                     beta =  beta,
                     E = E, #eV
                     Npoints=NX)
    z, x = kino.z()
    print('Focal length: {} m'.format(kino.f))
    
    # generate quadrant of radially symmetric pattern
    zq,xyq = radialSymmetric(z,x)
    idx = np.concatenate((-xyq[0],xyq[0]))
    del xyq
    print('Kinoform defined over square of size {} m \
          [{} to {} m, {} pixels]'.format( np.max(idx)-np.min(idx),
                                           np.min(idx),
                                           np.max(idx),                                         
                                           len(idx)))
    print('Thickness: {} m'.format(np.max(z)))
    
    
    # assemble 4 quadrants from the first one
    R=np.concatenate((np.rot90(zq),zq),0)
    C= np.concatenate((np.fliplr(R),R),1)
    
    if display == True:
        plt.imshow(C)
    
    if savePath is not None:
        np.save(savePath, C)
        print('Wrote {}'.format(savePath))
    
    #imageio.imwrite('kinoform2D.tif',np.uint8(255*C/np.max(C)))
       
    return C,idx
    

    
def distanceFromMedian2D(x):
    # return 2D square array of distances from centre of a mesh constructed from
    # 1D list of values x.  The centre is taken as the [m,m] where m = median(x)
   
    from scipy.spatial import distance
    n = len(x)
    origin = [[0,0]] #[[np.median(x),np.median(x)]]
    xy = [ [xi, yi] for xi in x for yi in x]
    dist = distance.cdist(xy,origin,"euclidean")
    return np.reshape(dist,(n,n))
    


def testMake2D(diameter=15.0e-6, width = 30e-6, NumberOfZones=36,pixelSize=50e-9):
    resolution = pixelSize
    pady, padx = 10, 100
    NX = int(diameter/resolution)
    NY=int(width/resolution) #int(NX/10)
    print('Nx, Ny ={},{}'.format(NX,NY))
    kino = kinoform(diameter=diameter, 
                     Nzones = NumberOfZones, 
                     delta = 0.00314069726,  # for Ni at 600 eV (0.0031406)
                     beta = 0.000591546006,   # for Ni at 1000 eV
                     E = 0.600, #eV
                     Npoints=NX)
    z, x = kino.z_symmetric()
    profile2D = make2D(z, (NY,1))
    profile2D = np.pad(profile2D, [(pady, ), (padx, )], mode='constant')
    plt.imshow(profile2D)
    
    np.save('profile2D.npy', profile2D)


    


def transmission(thickness,E=600, Element='Ni', density=None, rho=None, mu=None, verbose=False):
    # compute transmission through object with thickness (in m) for a given photon energy E in eV and element specified by atomic number Z
    # For simplicity, this only works for elements.  Can be expanded to compounds using xraydb
    
    thickness = thickness*100  # to convert from m to cm
     
    if rho is None:
        print("Rho lookup not yet supported, must specify value")
        return -1
        
    mat = rm.Material(Element,table='Henke',rho=rho)
    
    if mu is None:
       mu = [mat.get_absorption_coefficient(e) for e in E]
     
    T = [np.exp(-(m*thickness)) for m in mu]
    
    if verbose:
        print('Transmission={}, Thickness={}, Energy={}, Element={}, mu={}, rho={}'.format(T, thickness, E, Element, mu, rho))
    
    return T


def plotTransmission(thickness, E1, E2, N, Element):
    import matplotlib.pyplot as plt
    E = np.linspace(E1,E2,N)
    T = transmission(thickness,E=E,Element=Element)
    plt.plot(E, T)

def transmissionInt2D(thickness,E, Element='Ni', rho=8.902):
    
    # get the transmission in each thickness element
    tr = [transmission(t, E=E, Element=Element,rho=rho) for t in np.nditer(thickness) ]
    
    # get the average as the sum of transmission over all thickness elements divided by number of elements
    intTr = np.sum( tr)  / np.ndarray.size(thickness)
     
    return intTr

def test():
    
    fileThickness = '/user/home/wpg/extensions/out/check transmission/thickness-cropped.tif'
    thickness = imageio.imread(fileThickness)
    tr =  transmissionInt2D(thickness,E=600,Element='Ni')
    print (tr)
    
    
def phaseEff(E, E_0, Element, rho):
    mat = rm.Material(Element,table='Henke',rho=rho)
    delta_0 = 1-np.abs(mat.get_refractive_index(E_0)) 
    print(delta_0)
    delta = [1-np.real(mat.get_refractive_index(Ei)) for Ei in E]   
    eff =  [np.sinc(np.pi*(delta[i]/delta_0 * E[i]/E_0-1))**2    for i in range(len(E))]
    return eff




Element = 'Ni'
rho = 8.902  # g/cm^3

#Element = 'C'
#rho = 3.51 

#Element = 'Si'
#rho = 2.32

#Element='Si3N4'
#rho = 3.17

#Element ='Au'
#rho = 19.3

diameter=15.0e-6
NumberOfZones=38
pixelSize=2e-9
resolution = pixelSize
NX = int(diameter/resolution)

#mat = rm.Material(Element,rho=rho, table='Henke')
#mu = mat.get_absorption_coefficient(E)
#delta = 1.00-np.abs(mat.get_refractive_index(E))
#beta = -np.angle(mat.get_refractive_index(E)) 
Edb = np.loadtxt('xrayNi.dat',skiprows=2)
E = np.linspace(200,800,601, endpoint=True)
deltaL=np.interp(E, Edb[:,0], Edb[:,1])
betaL=np.interp( E, Edb[:,0], Edb[:,2])
wl = (1240./E)/1.e9



def beta(Eval=None):
    if Eval is not None:
        return betaL[np.where(E==Eval)][0]
    else:
        return betaL
    

def delta(Eval=None):
    if Eval is not None:
        return deltaL[np.where(E==Eval)][0]
    else:
        return deltaL


designEnergy = [300, 400, 500, 600, 700]

summary=[]
data=[]
for E0 in designEnergy:
       
    
    print('+++++++++++++++++++++++++++++++ Design Energy = {}'.format(E0))
    
    # initialise kinoform object
    kino = kinoform(diameter=diameter,  
                    Nzones = NumberOfZones,  
                    delta = delta(E0),#1-np.abs(mat.get_refractive_index(E0+DO)) ,  # for Ni at 600 eV (0.0031406) 
                    beta = beta(E0), #-np.angle(mat.get_refractive_index(E0+DO)),   # for Ni at 1000 eV  
                    E = E0, #eV
                    Npoints=NX)
    
    # get the radial positions and corresponding thickness values zz
    zz, xx = kino.z() 
    thickness = np.multiply(zz,100000)  # convert thickness profile to units of cm, scale by S
    #thickness = np.divide(thickness, 2)
    #thickness = [3e-5 for i in thickness]
    #thickness=[300e-5] #thickness[1:10]
    
    Tr = []
    #for m in mu:
    for b, w in zip(beta(),wl):    
        T=[]
        for th in thickness:
              #T.append (np.exp(-(m*th)))          
              T.append(np.exp(-(4*np.pi*b/w)*(th/100)))    
        Tr.append(np.sum(T)/len(T))
        
     
    
    #idealEff =  phaseEff(E,E0+DO, Element,rho=rho)
    idealEff = []
    for  Ei in E:
            #print ('delta0={}, delta={}, E={}'.format(delta(E0), delta(Ei), Ei ))
            arg = np.pi*( delta(Ei)/delta(E0) * Ei/(E0+0.0000000000001) -1)
            idealEff.append( np.square(np.sin(arg)/arg))     #idealEff = (np.sinc(np.pi*(( delta()/delta(E0)) * (E/E0) -1)))**2
    Efficiency = np.multiply(idealEff,Tr)
    
    # print summary information
    maxIdealEfficiency = E[np.where(idealEff==np.max(idealEff))]
    maxEfficiency = E[np.where(Efficiency == np.max(Efficiency))]
    
    print('Design energy = {} eV, \
          Design Thickness = {} nm'.format(E0,
                                           np.max(thickness)*10000000 ))
    print('Peak in efficiency {}: with absorption = {} eV, \n \
          without absorption = {} eV, \n \
          Shift = {} eV'.format(np.max(Efficiency), maxEfficiency,
                                maxIdealEfficiency,
                                maxEfficiency-maxIdealEfficiency))
    
    summary.append([E0,
                    np.max(thickness)*10000000,
                    np.max(Efficiency),
                    maxEfficiency,
                    np.max(idealEff),
                    maxIdealEfficiency,
                    maxIdealEfficiency-maxEfficiency])
    data.append([Tr, idealEff, Efficiency])
    
    # just plotting below
    fig, ax = plt.subplots(4,figsize=(5,15))
    fig.suptitle('Design energy = {} eV'.format(E0), fontsize=16)
    ax[0].plot(E,Tr,label='Transmission')
    ax[0].plot(E,idealEff,label='Eff. without absorption')
    ax[0].plot(E,Efficiency, label='Eff. with absorption')   
    #ax[0].xlim(np.min(E),np.max(E))
    ax[0].set(xlabel='Energy/eV')
    ax[0].legend()
    
    #ax[1].plot(E,np.divide(delta,beta),label='delta/beta')
    #ax[1].set(xlabel='Energy/eV')
    #ax[1].legend()
    
    ax[1].plot(E,Tr,label='Transmission')
    ax[1].set(xlabel='Energy/eV')
    #ax[1].set(yscale='log')
    ax[1].legend()
    
    
    ax[2].plot(E, delta(),label='delta')
    ax[2].plot(E, beta(), label='beta')
    ax[2].set(xlabel='Energy/eV')
    ax[2].set(yscale='log')
    ax[2].legend()
    
    ax[3].plot(np.multiply(xx,1000000),np.multiply(thickness,10000000))
    ax[3].set(xlabel='Radial distance [um]')
    ax[3].set(ylabel='Thickness [nm] ')
    ax[3].legend()  
    plt.show()
    
    
np.savetxt('summary.txt',summary)   

fig2, ax2 = plt.subplots(2,figsize=(10,10))
#ax2[0].plot(E,data[i][1], label='Efficiency without absorption')
for i in range(len(data)):
    ax2[0].plot(E,data[i][2],label='{} eV, thickness={} nm'.format(summary[i][1],summary[i][2]))
ax2[0].legend()
ax2[0].set_title('Efficiency with Absorption for different offset form design energy')
ax2[0].set(xlabel='Energy/eV')

#ax2[1].plot(np.min(E),np.min(data[0][0])) 
for i in range(len(data)):    
    ax2[1].plot(E,data[i][0],label='{} eV, thickness={} nm'.format(summary[i][1],summary[i][2]))
ax2[1].legend()
ax2[1].set_title('Transmission for different offset form design energy')
ax2[1].set(xlabel='Energy/eV')
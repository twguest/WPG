#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 06:28:17 2019

@author: gvanriessen



Note that the diameter parameter is not appropriately named, or there is a factor of 2
missing somewhere!  Needs to be checked!


"""

import numpy as np
from numpy import argmin
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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
    ax.plot(kino.x, kino.zEquiv()/scale, label='Equivalent profile / {}'.format(kino.Npoints))
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
    
    
def kinoform2D(diameter=15.0e-6, NumberOfZones=36,pixelSize=50e-9,savePath=None, display=False):

    resolution = pixelSize
    NX = int(diameter/resolution)*2
   
    kino = kinoform(diameter=diameter, 
                     Nzones = NumberOfZones, 
                     delta = 0.00314069726,  # for Ni at 600 eV (0.0031406)
                     beta = 0.000591546006,   # for Ni at 1000 eV
                     E = 0.600, #eV
                     Npoints=NX)
    z, x = kino.z()
    # generate quadrant of radially symmetric pattern
    zq,xyq = radialSymmetric(z,x)
    del xyq
    
    # assemble 4 quadrants from the first one
    R=np.concatenate((np.rot90(zq),zq),0)
    C= np.concatenate((np.fliplr(R),R),1)
    
    if display == True:
        plt.imshow(C)
    
    if savePath is not None:
        np.save(savePath, C)
    
    #imageio.imwrite('kinoform2D.tif',np.uint8(255*C/np.max(C)))
       
    return C, x
    

    
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
    
    
if __name__ == "__main__":
    #test(diameter=15.0e-6, NumberOfZones=36,pixelSize=50e-9)
    #compare()
    #testMake2D(diameter=15.0e-6, width = 30e-6, NumberOfZones=36,pixelSize=50e-9)
    kino, x  =   kinoform2D(diameter=15.0e-6, NumberOfZones=36,pixelSize=10.e-9, savePath='kinoformR30umPx10nm.npy')
    
    
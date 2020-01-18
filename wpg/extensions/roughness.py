#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:33:59 2019

@author: gvanriessen
"""
import numpy as np
from numpy.fft import ifft2, fft2, irfft2, rfft2
import matplotlib.pyplot as plt

def rsgeng2D(N,rL, h,clx=None,cly=None):
# adapted from rseng2D function distributed by  David Bergström (david.bergstrom@mysimlabs.com)
# =============================================================================
#     
# generates a square 2-dimensional random rough surface f(x,y) with NxN 
# surface points. The surface has a Gaussian height distribution function 
# and Gaussian autocovariance functions (in both x and y), where rL is the 
# length of the surface side, h is the RMS height and clx and cly are the 
# correlation lengths in x and y. Omitting cly makes the surface isotropic.
# 
# Input:    N   - number of surface points (along square side)
#           rL  - length of surface (along square side)
#           h   - rms height
#           clx, (cly)  - correlation lengths (in x and y)
# Output:   f  - surface heights
#           x  - surface points
#            y  - surface points
# =============================================================================

# =============================================================================
#      Nyquist sampling theorem sets a lower limit of the sample frequency of 
#       the surface to achieve the correct correlation length, see above (in 
#       this case the sampling frequency is referred to the sampling along 
#       one of the square sides of the surface). To achieve proper gaussian 
#       statistics in the 2D case use a ratio of surface length to correlation 
#       length rL/cl > ~70. Use the analysis code to check the quality of your 
#       surface and verify that you got what you asked for.
# 
# =============================================================================
    
## Note anisotropic surfaces not tested - note small imagineary part from ifft
     
    x = np.linspace(-rL/2,rL,N); 
    y = np.linspace(-rL/2,rL,N);
    [X,Y] = np.meshgrid(x,y); 

 
    Z = h*np.random.randn(N,N); # uncorrelated Gaussian random rough surface distribution
                      # with mean 0 and standard deviation h

    if clx is None and cly is None:
        f = Z

    # isotropic surface
    if clx is not None and cly is None:
    
        # Gaussian filter
        F = np.exp(-((X**2+Y**2)/(clx**2/2)));
    
        # correlation of surface including convolution (faltung), inverse
        # Fourier transform and normalizing prefactors
        #f = 2/sqrt(pi)*rL/N/clx*ifft2(fft2(Z).*fft2(F));
        f = s/np.sqrt(np.pi)*rL/N/clx * np.multiply( irfft2(fft2(Z),fft2(F)))
        
    # non-isotropic surface 
    if clx is not None and cly is not None:
    
        # Gaussian filter
        F = np.exp(-(X**2/(clx**2/2)+Y**2/(cly**2/2)));
    
        # correlated surface generation including convolution (faltning) and inverse
        # Fourier transform and normalizing prefactors
        A = np.real_if_close(ifft2(fft2(Z)),tol=1)
    
        B = np.real_if_close(fft2(F),tol=1)
        M = np.multiply(A,B)
        f = 2/np.sqrt(np.pi)*rL/N/np.sqrt(clx)/np.sqrt(cly)*M
        
    
    rmsf = rms(f)  # rms roughness
    
    return f, rmsf


def rms(x, axis=None):
    # return RMS of x, where x is a square array
    return sqrt(mean(x**2, axis=axis))

def testRoughness2D():
    lateralRes = 1.e-9   # should be small compared to roughness coherence lengths clx, cly
    h = 5e-9  
    rL = 30e-6  # ZP diameter
    N=rL/lateralRes
    rsgeng2d(N,rL,h)
    
    
    
    
    
    
    
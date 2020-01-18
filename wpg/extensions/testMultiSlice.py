#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:42:04 2019

@author: gvanriessen
"""

from wpg.optical_elements import _save_object, _load_object, Empty
from wpg.srwlib import srwl_uti_ph_en_conv
from wpg.generators import build_gauss_wavefront_xy
from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, plot_intensity_map, averaged_intensity,calculate_fwhm, get_intensity_on_axis
from wpg.useful_code.wfrutils import propagate_wavefront, calculate_fwhm_y,calculate_fwhm_x, print_beamline, get_mesh
from extensions.gvrutils import calculate_theta_fwhm_cdr
from extensions.multisliceOptE import *
import imageio


def constructTestWaveField(npoints=256,ekeV=0.6,z1=5, show=False):  #Alale: npoints = 256, z was 5
    # Generate a simple Gaussian wavefield for testing purposes

    # set wavefield parameters
    qnC        = 0.6
    wlambda    = srwl_uti_ph_en_conv(ekeV, _in_u='keV', _out_u='nm')
    theta_fwhm = calculate_theta_fwhm_cdr(ekeV,qnC)
    k          = 2*np.sqrt(2*np.log(2))
    range_xy   = (theta_fwhm/k*z1*5.)*2.0
    sigX       = 12.4e-10*k/(ekeV*4*np.pi*theta_fwhm)

    # construct wavefield
    wf0=build_gauss_wavefront_xy(nx=npoints, ny=npoints, ekev=ekeV,
                                              xMin=-range_xy/2 ,xMax=range_xy/2,
                                              yMin=-range_xy/2, yMax=range_xy/2,
                                              sigX=sigX, sigY=sigX,
                                              d2waist=z1,
                                              _mx=0, _my=0 )
    wfr = Wavefront(srwl_wavefront=wf0)    
    print('Wavelength=%f, theta FWWM=%f, range XY = %f, sig X = %f' % (wlambda, theta_fwhm, range_xy, sigX))

    wfr._srwl_wf.unitElFld = 1#'sqrt(Phot/s/0.1%bw/mm^2)'

    if show==True: #display the wavefield
        plotWavefront(wfr, 'Wavefield at source')

    return wfr


wf = constructTestWaveField(npoints=128)

print('Testing greyscale Slicer')
N=5
apple= greyscaleToSlices('extensions/in/blob.tif',slices=N,invert=False)
    
#apple.tr can be inserted into a beamline
apple.tr.propagate(wf)



print('Testing multisliceOptE')





wf = constructTestWaveField
slices=testblObjectSliced(6,trFilePath='/opt/wpg/wpg/extensions/in/Kinoform2978P.tif')#bZP30um_10nmres_200outerMostZoneWidth.tif')#30umLTZPAngleZero10nmRes.tif')#30um63D10nmRes.tif')#30umLTZP10nmRes63dUsingdataAtZeroAngle51pInXC.tif') #'./in/')  #30umLTZP10nmRes63dUsingdataAtZeroAngle51pInX.tif
multiSliceLTZP_X = multisliceOptE(objectSlices=slices,description='LTZP (multislice)')
multiSliceLTZP_X.printSliceParams()
multiSliceLTZP_X.description = 'LTZP X'
multiSliceLTZP_X.propagate(wf)
#cplx = multiSliceLTZP_X.getTransmissionFunction()
#phase = multiSliceLTZP_X.getTransmissionFunction(part = 'phase')
#amp = multiSliceLTZP_X.getTransmissionFunction(part='amplitude')

#multiSliceLTZP_X.showTransmissionFunction()

# now save tr
#multiSliceLTZP_X.writeTransmissionFunction('outx.pickle', 
#                                           comments='this is my transmission function, picled SRWLOPTT type') 

energy = 0.6 #keV
wavelength    = srwl_uti_ph_en_conv(energy, _in_u='keV', _out_u='nm') *1.0e-9 #m
transmission = multiSliceLTZP_X.getAmplitudeTransmission()
thickness  = multiSliceLTZP_X.getThickness()
phaseShift = multiSliceLTZP_X.getPhaseShift(wavelength)


import matplotlib.pyplot as plt
cmap = plt.get_cmap('PiYG')


fig, (ax0, ax1,ax2) = plt.subplots(nrows=3)

th = ax0.imshow(thickness,cmap=cmap)
fig.colorbar(th, ax=ax0)
ax0.set_title('Thickness')

ps = ax1.imshow(phaseShift,          cmap=cmap)
fig.colorbar(ps, ax=ax1)
ax1.set_title('Phase Shift')

tm = ax2.imshow(transmission,          cmap=cmap)
fig.colorbar(tm, ax=ax2)
ax2.set_title('Transmission')


fig.tight_layout()

plt.show()



imageio.imwrite('extensions/out/phaseShift.tif', np.float32(phaseShift))
imageio.imwrite('extensions/out/thickness.tif', np.float32(thickness))
imageio.imwrite('extensions/out/transmission.tif', np.float32(transmission))



# Test Plot - not ready for use
#from utilPlot import * 
#plotTr(multSliceLTZP_X.trOpt, mmin=0, mmax=1, pmin=-pi, pmax=+pi)

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

import matplotlib.pyplot as plt 




def constructTestWaveField(npoints=256,ekeV=0.6,z1=15, show=False):  #Alale: npoints = 256, z was 5
    # Generate a simple Gaussian wavefield for testing purposes

    # set wavefield parameters
    qnC        = 0.01
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



def testblObjectSliced(N,trFilePath,
                       resolution = [10.0e-9,10.0e-9,10.0e-9],
                       thickness= 5.e-7,       
                       delta= 4.7e-2,           #delta value of material in slice
                       attenLength= 0.277984e-5   #attenuation length in material in slice
                       ):
    """
     Function testing/demonstrating the generation of a list of N slices
     in structure expected for multisliceOptE

    """
    rx,ry,rz = resolution
    material=[thickness,delta,attenLength]


    zStep = 356.4e-9 
    lateralStep = 181.6e-9
 
    pp = propagationParameters(SemiAnalyt=1.0)

    slices=[]
    for i in range(N):
        zi = i* zStep  
        yi = 0.0
        xi = i* lateralStep -(lateralStep*N/2.0)
        r = 0.0
        pad =  None
        trans = [xi,yi,zi,r]
        slices.append( [trFilePath,trans,resolution,material,pad,pp] )

    return slices

def generate(wavefront, sourceImage = '/opt/wpg/wpg/extensions/in/blob.tif'):
    
    slices=testblObjectSliced(4,trFilePath=sourceImage)
    multiSlice = multisliceOptE(objectSlices=slices,description='Test (multislice)')
    multiSlice.printSliceParams()
    multiSlice.description = 'Test '
    multiSlice.propagate(wavefront)
    
    
    return multiSlice, wavefront



def getResults(tr):
    
    energy = 0.6 #keV
    wavelength    = srwl_uti_ph_en_conv(energy, _in_u='keV', _out_u='nm') *1.0e-9 #m
    
    transmission = tr.getAmplitudeTransmission()
    thickness  = tr.getThickness()
    phaseShift = tr.getPhaseShift(wavelength)
  
    return transmission, thickness, phaseShift



def plot(transmission, thickness, phaseshift):
    cmap = plt.get_cmap('PiYG')

    fig, (ax0, ax1,ax2) = plt.subplots(nrows=3)    

    th = ax0.imshow(thickness,cmap=cmap)
    fig.colorbar(th, ax=ax0)
    ax0.set_title('Thickness')
    
    ps = ax1.imshow(phaseshift,          cmap=cmap)
    fig.colorbar(ps, ax=ax1)
    ax1.set_title('Phase Shift')
    
    tm = ax2.imshow(transmission,          cmap=cmap)
    fig.colorbar(tm, ax=ax2)
    ax2.set_title('Transmission')  
    
    fig.tight_layout()
    plt.show()


def save(transmission, thickness, phaseShift):
    imageio.imwrite('out/phaseShiftTest.tif', np.float32(phaseShift))
    imageio.imwrite('out/thicknessTest.tif', np.float32(thickness))
    imageio.imwrite('out/transmissionTest.tif', np.float32(transmission))


def plotTr(tr):
    # Don't use this!  Trnsmission function is not comprised of a+bi type values
    from extensions.utilPlot import plotTr
    plotTr(tr.getSRWTransmissionFunction())#, mmin=0, mmax=1, pmin=-pi, pmax=+pi)

def plotThickness(tr):
    from extensions.utilPlot import plotTrThickness
    plotTrThickness(tr)
    
    
def plotPhaseShift(tr,wavelength):
    from extensions.utilPlot import plotTrPhaseShift
    plotTrPhaseShift(tr, wavelength)
    
'''class EpiCycleScalarFormatter(ticker.ScalarFormatter):
    def _set_orderOfMagnitude(self, rng):
        # can just be super() in py3, args only needed in LPy
        super(EpiCycleScalarFormatter, self)._set_orderOfMagnitude(rng)
        if self.orderOfMagnitude != 0:
            self.orderOfMagnitude -= 1
'''

    
    
#def test():
    
if 1>0:    
    ekeV = 600
    wlambda    = srwl_uti_ph_en_conv(ekeV, _in_u='keV', _out_u='nm')
    wf0 = constructTestWaveField(npoints=256,ekeV=ekeV,show=True)

    tr, wf1 = generate(wf0)
    transmission, thickness, phaseShift  = getResults(tr)
    plot(transmission, thickness, phaseShift)
    
          
    plotThickness(tr)
    
    from extensions.utilPlot import plotWavefieldComplex, multiple_formatter, Multiple
    from extensions.utilPlot import tickformat2, cb
    from matplotlib.ticker import FuncFormatter
    import cplot



    wft= constructTestWaveField(npoints=256,ekeV=ekeV,show=True)
    plotWavefieldComplex(wft)
    
    
'''
    #Problems with following - NAN valuess in phasesshift
    #plotPhaseShift(tr,wlambda)   
    
    wf0 = constructTestWaveField(npoints=256,ekeV=ekeV,show=True)


    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 14}

    plt.rc('font', **font)


    wf=wf0
    mesh = wf.params.Mesh
    nx, ny = mesh.nx,  mesh.ny
    xmin, xmax = mesh.xMin, mesh.xMax
    ymin, ymax = mesh.yMin, mesh.yMax

    b = wf.get_imag_part(slice_number=0)#[:,:,0]
    a = wf.get_real_part(slice_number=0)#[:,:,0]
    val = ne.evaluate("complex(a,b)")
    
    hx = (xmax - xmin) / nx
    x = np.linspace(xmin + hx / 2, xmax - hx / 2, nx)
    hy = (ymax - ymin) / ny
    y = np.linspace(ymin + hy / 2, ymax - hy / 2, ny)

    angle = np.arctan2(val.imag, val.real)
    #absval_scaled = np.abs(val)/np.max(np.abs(val)) #abs_scaling(np.abs(val))
    #srgb_vals = cplot.main.get_srgb(angle, absval_scaled)
    absval_scaled = np.abs(a)/np.max(np.abs(a))
    srgb_vals = cplot.main.get_srgb(b, absval_scaled)

    #setup axes with gridspec
    fig = plt.figure(figsize=(8,6))
    grid = plt.GridSpec(4,4, hspace=0.3, wspace=0.3)
    main_ax = fig.add_subplot(grid[0:2,1:3])
    #y_prof = fig.add_subplot(grid[0:2,0:1],xticklabels=[],sharey=main_ax)
    #x_prof = fig.add_subplot(grid[2,1:3], yticklabels=[], sharex=main_ax)
    cbf = fig.add_subplot(grid[0:2,3])#, yticklabels=[], xticklabels=[])
    cbf.set_xlabel('mag.')
    cbf.set_ylabel('phase (rad.)')
    
    im_main = main_ax.imshow( srgb_vals,
                              extent=(x.min(), x.max(), y.max(), y.min()),
                              interpolation="none",
                              origin="lower",
                              aspect="equal"
                             )
    
    main_ax.tick_params(axis='both', which='major')
    fmtx = FuncFormatter(lambda x, pos: tickformat2(x * 1.e6))
    fmty = FuncFormatter(lambda y, pos: tickformat2(y * 1.e6))
    main_ax.xaxis.set_major_formatter(fmtx)
    main_ax.yaxis.set_major_formatter(fmty)
    
    main_ax.set_xlabel('x [$\mu$m]')
    main_ax.set_ylabel('y [$\mu$m]')
   
    # now prepare and format colorar
    mmin,mmax=np.min(absval_scaled),np.max(absval_scaled)
    pmin,pmax=np.min(angle),np.max(angle)

    print('Phase: [%3.3f,%3.3f]' % (pmin,pmax))
    print('Abs: [%3.3f,%3.3f]' % (mmin,mmax))
    
    nc=200  # number of color levels
    cval = cb(mmin, mmax,pmin,pmax,nc)
        
    hx = (pmax - pmin) / nc
    x = np.linspace(pmin + hx / 2, pmax - hx / 2, nc)
    hy = (mmax - mmin) / nc
    y = np.linspace(mmin + hy / 2, mmax - hy / 2, nc)  
    
    angle = np.arctan2(cval.imag, cval.real)
    #absval_scaled = np.abs(cval) / np.max(np.abs(cval))
    srgb_cvals = cplot.main.get_srgb(cval.imag, cval.real)
    
    cbf.imshow(
        srgb_cvals,
        extent=(mmin, mmax, pmin, pmax),
        interpolation="none",
        origin="lower",
        aspect=0.75
    )
    
    cbf.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    cbf.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    cbf.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
    #cbf.grid(color='w', linestyle='-',linewidth=1)

    # plt.show()

'''

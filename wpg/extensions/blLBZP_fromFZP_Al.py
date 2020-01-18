#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:31:10 2018

@author: gvanriessen



In spyder run with:
    
    runfile('/opt/wpg/wpg/extensions/blLTZPeff4.1_Al.py', wdir='/opt/wpg/wpg/extensions', args='BLParams')
    
    or
    
    runfile('/opt/wpg/wpg/extensions/blLTZPeff4.1_Al.py', wdir='/opt/wpg/wpg/', args='extensions.BLParams')
"""

from wpg.optical_elements import _save_object, _load_object, Empty
from wpg.srwlib import srwl_uti_ph_en_conv
from wpg.generators import build_gauss_wavefront_xy
from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, plot_intensity_map, averaged_intensity,calculate_fwhm, get_intensity_on_axis
from wpg.useful_code.wfrutils import propagate_wavefront, calculate_fwhm_y,calculate_fwhm_x, print_beamline, get_mesh

import os, sys, imageio
from slugify import slugify
import matplotlib.pyplot as plt  # required?
import importlib
import time

from extensions.gvrutils import propagationParameters, calculate_theta_fwhm_cdr, \
                                plotWavefront,plot_wfront2
from extensions.gvrutils import writeIntensity, pixelScale
from extensions.gvrutils import get_transmissionFunctionCplx, calculate_fwhm
from extensions.gvrutils import get_intensity_on_axis
from extensions.gvrutils import show_transmissionFunctionCplx,  visVol,Struct
from extensions.multisliceOptE import *  

from extensions.beamlineExt import *


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

def testLTZPBeamline(conf):
# =============================================================================
#      test propgation through a  model beamline containing LTZP
#
# =============================================================================
    import datetime, uuid
  
    p = dynamic_import(conf)
    
    # make a unique output path
    strOutPath = p.OUTPUTPATHBASE + datetime.date.today().isoformat() +'/'+  str(uuid.uuid1()) + '/'
    try: 
        os.makedirs(strOutPath)
        print('Output directory created [%s]' % strOutPath) 
    except OSError:
        SAVE = False
        print('Problem creating output directory, output will NOT BE SAVED!!!')
        pass

    volOutPath =  strOutPath+'/vol/'

    blSaveInitial   = oeScreen(strOutPath   + p.wf0Path,   
                                  description='Initial WF')
    blSavePreLTZPX  = oeScreen(strOutPath + p.wfPreOpticPath, 
                                  description='Before LTZP',
                                  writePhase=True,
                                  ROI = p.preLTZPROI ) 
    blSavePostLTZPX = oeScreen(strOutPath + p.wfPostZPXPath,   
                                  description='After ZPX',
                                  showPhase=False,
                                  writePhase=True)
    blSaveFocus     = oeScreen(strOutPath     + p.wfFocusXPath,
                                  description='At Focus',
                                  showPhase=False,
                                  writePhase=True)
    blSaveDetector  = oeScreen(strOutPath  + p.wfDetectorPath, 
                             description='At detector')
    blSavePostOSA  = oeScreen(strOutPath     + p.wfPostOSAPath,
                                  description='Post OSA',
                                  ROI=p.postOSAROI,
                                  showPhase=False,
                                  writePhase=True)
#    blSavePostCS  = oeScreen(strOutPath     + p.wfPostCSPath,
#                                  description='Post CS',
#                                  ROI=p.postOSAROI,
#                                  showPhase=False,
#                                  writePhase=True)
#    
    # White beam slits  After the source 
    elWBS = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "a", p.WBS.xwidth, p.WBS.ywidth, 0.0, 0.0)], 
                        _arProp=[p.WBS.pp]),
               description = p.WBS.description)
   
    #drift 6.5m
    elDriftA = bl(SRWLOptC(_arOpt=[SRWLOptD(p.DriftA.driftLength)], 
                  _arProp=[p.DriftA.pp]), 
                  description= p.DriftA.description)

    # Beamline Slits.  These define a secondary source
    elBLS = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "a", p.BLS.xwidth, p.BLS.ywidth, 0.0, 0.0)],
               _arProp=[p.BLS.pp]),
               description=p.BLS.description)
               
    #drift 6.5m
    elDriftCleanup = bl(SRWLOptC(_arOpt=[SRWLOptD(p.DriftCleanup.driftLength)], 
                        _arProp=[p.DriftCleanup.pp]), 
                        description=p.DriftCleanup.description)     
    
    # Beamline Slits.  These define a secondary source
    elBLSCleanup = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "a", p.BLSCleanup.xwidth, p.BLSCleanup.ywidth, 0.0, 0.0)],
                      _arProp=[p.BLSCleanup.pp]),
                      description=p.BLSCleanup.description)
   
    #drift 6.5m
    elDriftB = bl(SRWLOptC(_arOpt=[SRWLOptD(p.DriftB.driftLength)], 
                  _arProp=[p.DriftB.pp]), 
                  p.DriftB.description)

    # Beam defining aperture
    elBDA = bl(SRWLOptC(_arOpt=[SRWLOptA("c", "a", p.BDA.xwidth, p.BDA.ywidth, 0.0, 0.0)],
               _arProp=[p.BDA.pp]),
               p.BDA.description)

    elDriftBDA = bl(SRWLOptC(_arOpt=[SRWLOptD(p.driftBDA.driftLength)], 
                    _arProp=[p.driftBDA.pp]), 
                    description=p.driftBDA.description)
    
    elResc = bl(SRWLOptC(_arOpt=[SRWLOptD(p.Resc.driftLength)], 
                    _arProp=[p.Resc.pp]), 
                    description=p.Resc.description)
    
    
    elCS = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "o", p.CS.xwidth, p.CS.ywidth, 0.0, 0.0)],
               _arProp=[p.CS.pp]),
               description=p.CS.description)
    

    blDriftToFocus =   bl(SRWLOptC(_arOpt=[SRWLOptD(p.DriftToFocus.driftLength)], 
                          _arProp=[p.DriftToFocus.pp]), 
                          description=p.DriftToFocus.description)
                        
    #from  extensions.driftVol import driftVol 
    #blfocusVol=driftVol(p.focusVol.distances, p.focusVol.pp, volOutPath) 
    #blfocusVol.description=p.focusVol.description 
    
    blRescAfterF = bl(SRWLOptC(_arOpt=[SRWLOptD(p.RescAfterF.driftLength)], 
                    _arProp=[p.RescAfterF.pp]), 
                    description=p.RescAfterF.description)
    
    # Order sorting aperture
    elOSA = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "a", p.OSA.xwidth, p.OSA.ywidth, 0.0, 0.0)],
               _arProp=[p.OSA.pp]),
               description=p.OSA.description)
    
    blRescAfterOSA = bl(SRWLOptC(_arOpt=[SRWLOptD(p.RescAfterOSA.driftLength)], 
                    _arProp=[p.RescAfterOSA.pp]), 
                    description=p.RescAfterOSA.description)

    
    zeroDrift = bl(SRWLOptC(_arOpt=[SRWLOptD(p.zeroDrift.driftLength)], 
                         _arProp=[p.zeroDrift.pp]),
                         description =  p.zeroDrift.description)
    

    driftToDetector = bl(SRWLOptC(_arOpt=[SRWLOptD(p.driftToDetector.driftLength)], 
                         _arProp=[p.driftToDetector.pp]),
                         description =  p.driftToDetector.description)
    #driftToDetector.absParam = ([driftToDetector.range,driftToDetector.range,driftToDetector.resol ,driftToDetector.resol ])

    
    if not( os.path.exists(p.CACHEPATH+p.wfPostZPXPath) or (os.path.exists(p.CACHEPATH+p.wfFocusXPath))):    
        # construct beamlines containing  test object oriented in x direction     
        
        blCleanupOptic = bl(SRWLOptC(_arOpt=[SRWLOptA("r", "a", p.LTZPx.aperture.xwidth, p.LTZPx.aperture.ywidth, 0.0, 0.0)],
                            _arProp=[p.LTZPx.aperture.pp]),
                            description=p.LTZPx.aperture.description)

       
        p.LTZPx.slices = testblObjectSliced(p.LTZPx.NSlices,
                                          p.LTZPx.trFile,
                                          orientation = p.LTZPx.orientation,
                                          resolution  = p.LTZPx.resolution,
                                          thickness   = p.LTZPx.thickness,
                                          delta       = p.LTZPx.delta,
                                          attenLength = p.LTZPx.attenLength,
                                          trPreFilePath = p.LTZPx.trPreFile  
                                         )
        
        multiSliceLTZP_X = multisliceOptE(objectSlices=p.LTZPx.slices,description=p.LTZPx.description)
        multiSliceLTZP_X.printSliceParams()
       
            
    if (os.path.exists(p.CACHEPATH+p.wfFocusXPath)):
        wfr = Wavefront()
        wfr.load_hdf5(p.CACHEPATH+p.wfFocusXPath)
        
        OpticalSystem = [blRescAfterF,
                         elOSA,
                       driftToDetector]
                      #blRescAfterOSA ]
                      #blSaveDetector]
    
    elif (os.path.exists(p.CACHEPATH+p.wfPostZPXPath)):
        wfr = Wavefront()
        wfr.load_hdf5(p.CACHEPATH+p.wfPostZPXPath)
        
        OpticalSystem = [
                     blDriftToFocus, #blSaveFocus
                     elOSA,
                     #blSavePostOSA,
                     driftToDetector]
                     #blSaveDetector]
                     
#        OpticalSystem = [                     
#                     focusVol,
#                     blSaveFocus,
#                     elOSA,
#                     blRescAfterOSA,
#                     driftToDetector,
#                     blSaveDetector]
                     
    elif (os.path.exists(p.CACHEPATH+p.wfPreOpticPath)):
        wfr = Wavefront()
        wfr.load_hdf5(p.CACHEPATH+p.wfPreOpticPath)
        
        OpticalSystem = [ #blCleanupOptic, #
                         elCS,
                        # blSavePostCS,                               
                         multiSliceLTZP_X,
                         blSavePostLTZPX,                        
                         blDriftToFocus,
                         blSaveFocus,
                         elOSA,
                         driftToDetector]
                        # blSaveDetector]

#        OpticalSystem = [
#                         blCleanupOptic,        
#                         multiSliceLTZP_X, 
#                         blSavePostLTZPX,                        
#                         focusVol,
#                          blSaveFocus,
#                         elOSA,
#                         blRescAfterOSA,
#                         driftToDetector,
#                         blSaveDetector]

    else:         
        # generate a generic Gaussian beam for testing
        wfr = constructTestWaveField(npoints=500, ekeV=0.6, show = p.SHOWPLOTS)
             
    #drift region between optics
       # blDriftXY = blDrift(50.e-6)
       # blDriftXY.description = 'Drift LTZPX to LTZPY'
        
    #        construct beamlines containing  test object oreinted in y direction
    #        slices = testblObjectSliced(N,trFile,
    #                                    orientation='x',
    #                                    resolution = [10.0e-9,10.0e-9,10.0e-9],
    #                                    thickness= 657.85e-9,
    #                                    delta= 0.0031406,
    #                                    attenLength= 0.277984e-6
    #                                    )
    
        #multiSliceLTZP_Y = multisliceOptE(objectSlices=slices,description='LTZP (multislice)')
        #multiSliceLTZP_Y.printSliceParams()
        #multiSliceLTZP_X.description = 'LTZP Y'
      
        OpticalSystem = [elWBS, 
                         elDriftA, 
                         elBLS, 
                         elDriftCleanup, 
                         elBLSCleanup, 
                         elDriftB, 
                         elBDA, 
                         elDriftBDA, 
                         elResc,
                         #blCleanupOptic, 
                         elCS,
                         blSavePreLTZPX,                         
                         multiSliceLTZP_X, 
                         blSavePostLTZPX,
                         blDriftToFocus,   # drftVol
                         blSaveFocus, 
                         elOSA, 
                         blSavePostOSA,
                         driftToDetector,
                         blSaveDetector]

    # Now propagate wavefield through all elemens in OpticalSysstem
    propResults = []
    for section in OpticalSystem:
                  
          # display summary of beamline  section
          #section.print_beamline()
          print ('Progagating through %s' % section.description)
          starttime = time.time()
          metrics = None
          try:
              metrics = section.propagate(wfr=wfr,  
                                            plot=p.SHOWPLOTS, 
                                            outFilePath=None, 
                                            propagationOverride=section.absParam,
                                            describe=False)  
                  
          except AttributeError:              
              try:
                  section.propagate(wfr=wfr,  plot=p.SHOWPLOTS, outFilePath=None)  
              except TypeError:
                  section.propagate(wfr=wfr)        
            
          if p.SAVE_SNAPSHOTS is True: #Save intensity image
              writeIntensity(wfr,
                             section.description,
                             path=strOutPath,
                             polarization='horizontal',
                             imgType='tif')
#             
              plot_map_I(wfr, save=strOutPath+ section.description+'_profile.tif', range_x=35.0e-3, range_y=35.0e-3,im_aspect='equal') #range_x in mm
          t = time.time() - starttime
          print('Propagation through {} took {} seconds'.format(section.description, t))
          propResults.append([section.description, t, metrics])

    #np.savetxt(strOutPath+'/summary.csv',propResults,delimeter=',')
    
    if multiSliceLTZP_X in OpticalSystem: 
        transmission = multiSliceLTZP_X.getAmplitudeTransmission()
        thickness  = multiSliceLTZP_X.getThickness()
        phaseShift = multiSliceLTZP_X.getPhaseShift(p.wavelength)

        cmap = plt.get_cmap('PiYG')
        
        fig, (ax0, ax1,ax2) = plt.subplots(nrows=3)
       
        th = ax0.imshow(thickness,cmap=cmap)
        fig.colorbar(th, ax=ax0)
        ax0.set_title('Thickness')
        
        ps = ax1.imshow(phaseShift,  cmap=cmap)
        fig.colorbar(ps, ax=ax1)
        ax1.set_title('Phase Shift')
        
        tm = ax2.imshow(transmission, cmap=cmap)
        fig.colorbar(tm, ax=ax2)
        ax2.set_title('Amplitude Transmission')
              
        fig.tight_layout()      
        plt.show(block = False)

        imageio.imwrite(strOutPath+'/phaseShift.tif', np.float32(phaseShift))
        imageio.imwrite(strOutPath+'/thickness.tif', np.float32(thickness))
        imageio.imwrite(strOutPath+'/transmission.tif', np.float32(transmission))
        
    # now save tr
    multiSliceLTZP_X.writeTransmissionFunction(strOutPath+'/LTZP63d_transmission.pickle', 
                                           comments='this is my transmission function, pickled SRWLOPTT type') 
  
     
                                          

 
def dynamic_import(module):
 
    return importlib.import_module(module)
 
def main(conf):
     
#     try:
#         import_file(conf)
#     except:
#         print('Could not load specified parameter file')     
#     else:
#        
     #import_file(conf)  
    # dynamic_import(conf)
     testLTZPBeamline(conf)

if __name__ == '__main__':
    
    try:
        print(sys.argv[1])
        main(sys.argv[1])
    except IndexError as e:
        #main('extensions.parameters.BLParams_forKino')
        main('extensions.parameters.BLParams_forLBZP')
        

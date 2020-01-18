#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:31:10 2018

@author: gvanriessen



In spyder run with:
    
    runfile('/opt/wpg/wpg/extensions/blKINOFORMeff4.1_Al.py', wdir='/opt/wpg/wpg/extensions', args='BLParams')
    
    or
    
    runfile('/opt/wpg/wpg/extensions/blKINOFORMeff4.1_Al.py', wdir='/opt/wpg/wpg/', args='extensions.BLParams')
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
import csv
import numpy as np

from extensions.gvrutils import propagationParameters, calculate_theta_fwhm_cdr, \
                                plotWavefront,plot_wfront2
from extensions.gvrutils import writeIntensity, pixelScale
#from extensions.gvrutils import get_transmissionFunctionCplx,
from extensions.gvrutils import  calculate_fwhm
from extensions.gvrutils import get_intensity_on_axis
from extensions.gvrutils import show_transmissionFunctionCplx,  visVol,Struct
from extensions.multisliceOptE import *  
from extensions.beamlineExt import *
from extensions.wf_sampling import getDimensions,  modifyRangeWF,modifyResolutionWF, modifyPixelsWF

from extensions.wavefrontGenerator import gsmSource


def beam(E):
    
    nx, ny = 6096, 6096  
    dimensions = [-1000e-6, 1000e-6, -1500e-6, 1500e-6] # dimensions
    gsm = gsmSource( 1, 1,        # number of modes in x and y
                    100e-6,1e-6,  #rho_I, rho_mu for x modes
                    100e-6, 1e-6,  #rho_I, rho_mu for y modes
                    E, # keV
                    -10, #radius of curvature,
                    nx, ny,  dimensions) # dimensions)
   
    gsm.generateModes(show=False)  #called explicitly, not done during init because slow
    
    gsm.getDimensions()
    
    return gsm
    


def writeMetrics(path, metrics, elementDescription):
    
    fmet = open(path,'a')  
    cs = csv.writer(fmet)
    cs.writerow([elementDescription])
    if metrics is not None:
        try:
            for key, val in metrics.items.items():
                    cs.writerow([key, val])
        except:
            cs.writerow(str(metrics))
    
    fmet.close()


#def constructTestWaveField(npoints=4096,
#                           ekeV=0.6,
#                           waist=-10,
#                           range_xy=2500e-6,
#                           show=False,
#                           theta_fwhm=100e-6): 
#    # Generate a simple Gaussian wavefield for testing purposes
#
#    # set wavefield parameters
#    wlambda    = srwl_uti_ph_en_conv(ekeV, _in_u='keV', _out_u='nm')    
#    k          = 2*np.sqrt(2*np.log(2))
#    #range_xy   = (theta_fwhm/k*z1*5.)*2.0
#    sigX       = 12.4e-10*k/(ekeV*4*np.pi*theta_fwhm)
#
#    # construct wavefield
#    wf0=build_gauss_wavefront_xy(nx=npoints, ny=npoints, ekev=ekeV,
#                                              xMin=-range_xy/2 ,xMax=range_xy/2,
#                                              yMin=-range_xy/2, yMax=range_xy/2,
#                                              sigX=sigX, sigY=sigX,
#                                              d2waist=waist,
#                                              _mx=0, _my=0 )
#    wfr = Wavefront(srwl_wavefront=wf0)    
#    print('Wavelength=%f, theta FWWM=%f, range XY = %f, sig X = %f' % (wlambda, theta_fwhm, range_xy, sigX))
#    
#    wfr._srwl_wf.unitElFld = 1#'sqrt(Phot/s/0.1%bw/mm^2)'
#    
#    if show==True: #display the wavefield
#        plotWavefront(wfr, 'Wavefield at source')
#
#    return wfr

def testKINOFORMBeamline(conf,roughScale=None):
# =============================================================================
#      test propgation through a  model beamline containing KINOFORM
#
# =============================================================================
    import datetime, uuid
  
    p = dynamic_import(conf)
    
    # scale mask roughness, overriding value in param file
    if roughScale is not None:
        p.mask.S = roughScale
        print ('#*#*#*#*#*# Roughscale = {} (type={})'.format(roughScale,type(roughScale )))
        scaledMaskThickness = p.mask.S*p.mask.h
        
        r = imageio.imread(p.mask.file)
        meanMaskThickness = (np.mean(r)/255)*scaledMaskThickness
        del r  # in future we should keep it and pass the numpy array to the optTr function
        
        scaledKinoThickness = p.KINOFORMx.thickness-meanMaskThickness
        
        print ('Mask thickness:  max={}, mean ={}\n Optic thickness: max.={}'.format(scaledMaskThickness, meanMaskThickness, scaledKinoThickness) ) 
        
        p.LABEL = "{}-RMS{}".format(p.LABEL,roughScale)
        p.OUTPUTPATHBASE = p.BASEPATH+'out/' + p.LABEL + '/'
        p.focusVol.outPath = p.OUTPUTPATHBASE+'/focusVol/'

    
    # make a unique output path
    strOutPath = p.OUTPUTPATHBASE + datetime.date.today().isoformat() +'/'+  str(uuid.uuid1()) + '/'
    try: 
        os.makedirs(strOutPath)        
        print('Output directory created [%s]' % strOutPath) 
    except OSError:
        SAVE = False
        print('Problem creating output directory, output will NOT BE SAVED!!!')
        
    try: 
        os.makedirs(strOutPath+'multislice')        
        print('Output directory created [%s]' % (strOutPath+'multislice')) 
    except OSError:
        SAVE = False
        print('Problem creating multislice output directory')
        

    EkeV = p.EkeV

    volOutPath =  strOutPath+'/vol/'

    blSaveInitial   = oeScreen(strOutPath   + p.wf0Path,   
                                  description='Initial WF')
    blSavePreKINOFORMX  = oeScreen(strOutPath + p.wfPreOpticPath, 
                                  description='Before KINOFORM',
                                  writePhase=True,
                                  ROI = p.preKINOFORMROI ) 
    blSavePostKINOFORMX = oeScreen(strOutPath + p.wfPostZPXPath,   
                                  description='After ZPX',
                                  showPhase=False,
                                  writePhase=True)
    blSavePostMask = oeScreen(strOutPath + p.wfPostMaskPath,   
                                  description='After Mask',
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
    blSavePostCS  = oeScreen(strOutPath     + p.wfPostCSPath,
                                  description='Post CS',
                                  ROI=p.postOSAROI,
                                  showPhase=False,
                                  writePhase=True)
    blSaveBeam  = oeScreen(strOutPath     + p.wfBeamPath,
                                  description='Beam (pre CS)',                                  
                                  showPhase=False,
                                  writePhase=True)
    
    
    
    # White beam slits  After the source 
    elWBS = bl(SRWLOptC(_arOpt=[SRWLOptA("c", "a", p.WBS.xwidth, p.WBS.ywidth, 0.0, 0.0)], 
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
    
    
    elCS = bl(SRWLOptC(_arOpt=[SRWLOptA("c", "o", p.CS.diameter, 0.0, 0.0)],
               _arProp=[p.CS.pp]),
               description=p.CS.description)
    
    maskTr = srwl_opt_setup_transm_from_file( p.mask.file,p.mask.resolution,
                                              scaledMaskThickness, 
                                              p.mask.delta, 
                                              p.mask.attenLength,
                                              background_color=255,
                                              invert=True)
    
    blMask = bl(SRWLOptC(_arOpt=[maskTr],
                         _arProp=[p.mask.pp]),
                         description=p.mask.description)
    

    blDriftToFocus =   bl(SRWLOptC(_arOpt=[SRWLOptD(p.DriftToFocus.driftLength)], 
                          _arProp=[p.DriftToFocus.pp]), 
                          description=p.DriftToFocus.description)
                        
    from  extensions.driftVol import driftVol 
  
    # define driftVol
    distances=(p.focusVol.distances)
    blfocusVol=driftVol(distances,
                  p.focusVol.pp,
                  p.focusVol.outPath, 
                  zoomXOutput = 0.25,     # reduce size of output, keep all significant features
                  zoomYOutput = 0.25,     # reduce size of output, keep all significant features
                  resampleOutput=2,     # reduce size of output, proportional loss of resolution
                  method=1,
                  saveIntensity=True,
                  saveComplex=False)
    blfocusVol.description=p.focusVol.description 
    
  
    blRescAfterF = bl(SRWLOptC(_arOpt=[SRWLOptD(p.RescAfterF.driftLength)], 
                    _arProp=[p.RescAfterF.pp]), 
                    description=p.RescAfterF.description)
    
    # Order sorting aperture
    elOSA = bl(SRWLOptC(_arOpt=[SRWLOptA("c", "a", p.OSA.diameter, 0.0, 0.0)],
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

    
    if p.doFocusVolume is True:
        blWaist =  blfocusVol
    else:
        blWaist = blDriftToFocus
    
    cacheA = os.path.exists(p.CACHEPATH+p.wfPostZPXPath)
    cacheB = os.path.exists(p.CACHEPATH+p.wfFocusXPath)
    print('Checking for {}... Found: {}'.format(p.CACHEPATH+p.wfPostZPXPath,cacheA ))
    print('Checking for {}... Found: {}'.format(p.CACHEPATH+p.wfFocusXPath,cacheB ))
    
    if cacheA is False and cacheB is False:    
        # construct beamlines containing  test object oriented in x direction    
        
        
        blCleanupOptic = bl(SRWLOptC(_arOpt=[SRWLOptA("c", "a", p.KINOFORMx.aperture.xwidth, p.KINOFORMx.aperture.ywidth, 0.0, 0.0)],
                            _arProp=[p.KINOFORMx.aperture.pp]),
                            description=p.KINOFORMx.aperture.description)

       
        p.KINOFORMx.slices = testblObjectSliced(p.KINOFORMx.NSlices,
                                          p.KINOFORMx.trFile,
                                          orientation = p.KINOFORMx.orientation,
                                          resolution  = p.KINOFORMx.resolution,
                                          thickness   = scaledKinoThickness, # compensate for mask
                                          delta       = p.KINOFORMx.delta,
                                          attenLength = p.KINOFORMx.attenLength,
                                          trPreFilePath = None,
                                          trPreProperties = None
                                         )
        
        multiSliceKINOFORM_X = multisliceOptE(objectSlices=p.KINOFORMx.slices,description=p.KINOFORMx.description)
        multiSliceKINOFORM_X.printSliceParams()
       
            
    if cacheB is True:
        wfr = Wavefront()
        wfr.load_hdf5(p.CACHEPATH+p.wfFocusXPath)
        
        OpticalSystem = [elOSA, 
                         blSavePostOSA,
                         driftToDetector,
                         blSaveDetector]
    
    elif (cacheA is True):
        wfr = Wavefront()
        wfr.load_hdf5(p.CACHEPATH+p.wfPostZPXPath)
        
        OpticalSystem = [blMask,
                         blSavePostMask,
                         blWaist,
                         blSaveFocus, 
                         elOSA, 
                         blSavePostOSA,
                         driftToDetector,
                         blSaveDetector]
                     
#    elif (os.path.exists(p.CACHEPATH+p.wfBeamPath)):
#        wfr = Wavefront()
#        wfr.load_hdf5(p.CACHEPATH+p.wfBeamPath)
#        
#        OpticalSystem = [#blCleanupOptic, #
#                         elCS,
#                         blSavePostCS,                               
#                         multiSliceKINOFORM_X,
#                         blSavePostKINOFORMX,                        
#                         blWaist,
#                         blSaveFocus,
#                         elOSA,
#                         blSavePostOSA,
#                         driftToDetector,
#                         blSaveDetector]


    else:         
        # generate a generic Gaussian beam for testing
        source = beam(p.EkeV); # constructTestWaveField(npoints=2048, ekeV=p.EkeV, show = p.SHOWPLOTS)
        wfr = source.modes[0]
        srcBeamline = [#elWBS, 
                         elDriftA,
                          blSaveInitial
                         #elBLS, 
                         #elDriftCleanup, 
                         #elBLSCleanup, 
                         #elDriftB, 
                         #elBDA, 
                         #elDriftBDA, 
                         # elResc]
        ]
        
        for section in srcBeamline:        
#             metrics = section.propagate(wfr=wfr,  
#                                            plot=p.SHOWPLOTS, 
#                                            outFilePath=None, 
#                                            #propagationOverride=section.absParam,
#                                            describe=True)  
             section.propagate(wfr)
             #writeMetrics(strOutPath+'/'+p.LABEL+'-summary.csv', metrics, section.description) 
        
        nx, ny, dx, dy, rx, ry = getDimensions(wfr)     
        # Range and resolution that is required for next step:
        Rx, Ry = 100e-6, 100e-6    #suspect problem with asymmetric wfr - why?
        Dx, Dy = p.KINOFORMx.resolution[0], p.KINOFORMx.resolution[1]
        print ('*********************************** \n Wavefront resolution before optic will be set to {} m x {} m ********************************************'.format(Dx,Dy) )
        
        modifyRangeWF(wfr, Rx, Ry)
        modifyResolutionWF(wfr, Dx, Dy)
        #modifyPixelsWF(wfr, 22000,22000)
   
        nx, ny, dx, dy, rx, ry = getDimensions(wfr)
        print('After resize and rescale: Number={}, {}, Range={}, {}, Step={}, {}'.format(nx,ny,rx,ry,dx,dy))

        blSaveBeam.propagate(wfr=wfr)


        OpticalSystem = [ #blCleanupOptic, 
                         elCS,
                         blCleanupOptic,
                         blSavePreKINOFORMX,                         
                         multiSliceKINOFORM_X, 
                         blSavePostKINOFORMX,
                         blMask,
                         blSavePostMask,
                         blWaist,
                         blSaveFocus, 
                         elOSA, 
                         blSavePostOSA,
                         driftToDetector,
                         blSaveDetector]

    # Now propagate wavefield through all elements in OpticalSysstem
        
    propResults = []
    for section in OpticalSystem:
                  
          # display summary of beamline  section
          #section.print_beamline()
          print ('Progagating through %s' % section.description)
          starttime = time.time()
          metrics = None
          
          if type(section) is oeScreen:
              metrics = section.propagate(wfr=wfr)
#              
#              #make a symbolic link to WF in cache
#              if p.updateCache:
#                  path=section.filename + section.extension
#                  if os.path.exists(path):
#                      os.rename(path,path+'.replaced')
#                  
#                  os.symlink(path, p.CACHEPATH + os.path.basename(path))
#                  
              
          
          elif isinstance(section,multisliceOptE):
                section.propagate(wfr)#, outFilePath=strOutPath+'multislice/' )
              
                transmission = multiSliceKINOFORM_X.getAmplitudeTransmission()
                thickness  = multiSliceKINOFORM_X.getThickness()
                phaseShift = multiSliceKINOFORM_X.getPhaseShift(p.wavelength)
        
                imageio.imwrite(strOutPath+'/phaseShift.tif', np.float32(phaseShift))
                imageio.imwrite(strOutPath+'/thickness.tif', np.float32(thickness))
                imageio.imwrite(strOutPath+'/transmission.tif', np.float32(transmission))
                
                # now save tr
                #multiSliceKINOFORM_X.writeTransmissionFunction(strOutPath+'/KINOFORM63d_transmission.pickle', 
                #                                   comments='this is my transmission function, pickled SRWLOPTT type') 
          
          elif isinstance(section,driftVol):
              section.propagate(wfr,  plot=False,outFilePath=None)
              
              
              
          else:   
              
              #handle differences between propagate methods for other object types
              try:
                  metrics = section.propagate(wfr=wfr,  
                                                plot=p.SHOWPLOTS, 
                                                outFilePath=None, 
                                                #propagationOverride=section.absParam,
                                                describe=True)  
                      
              except:# (AttributeError, TypeError) as err1:              
                  try:
                      section.propagate(wfr=wfr,  plot=p.SHOWPLOTS, outFilePath=None)  
                  except:  # (TypeError, AttributeError) as err2:
                      section.propagate(wfr=wfr)        
              
                
                
            
          if p.SAVE_SNAPSHOTS is True: #Save intensity image
              writeIntensity(wfr,
                             section.description,
                             path=strOutPath,
                             polarization='horizontal',
                             imgType='tif')             
              #plot_map_I(wfr, save=strOutPath+section.description+'_profile.tif', range_x=35.0e-3, range_y=30.0e-3,im_aspect='equal') #range_x in mm
          t = time.time() - starttime
          print('Propagation through {} took {} seconds'.format(section.description, t))
          propResults.append([section.description, t, metrics])



          writeMetrics(strOutPath+'/'+p.LABEL+'-summary.csv', metrics, section.description) 
          
    #np.savetxt(strOutPath+'/'+p.LABEL+'-summary.csv',propResults,delimiter=',')

    print ('Kinoform focal length = {}.  Mask roughness scaling factor = {}'.format(p.F1,p.mask.S))
    print('Output: {}'.format(strOutPath))

 
def dynamic_import(module):
 
    return importlib.import_module(module)
 
def main(conf, RMSScalingFactor):
     
#     try:
#         import_file(conf)
#     except:
#         print('Could not load specified parameter file')     
#     else:
#        
     #import_file(conf)  
    # dynamic_import(conf)
     plt.ioff()
     testKINOFORMBeamline(conf, RMSScalingFactor)

if __name__ == '__main__':
    

        main(sys.argv[1], float(sys.argv[2]))
#        if len(sys.argv) == 3:
#            print('Parameters {}, R={}'.format(sys.argv[1],sys.argv[2]))
#            main(sys.argv[1], sys.argv[2])
#        else: #assume execution for testing only
#            print('No arguments supplid, assume testing mode...')
#            main ('extensions.parameters.BLParams_Kino',1)
#        
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:53:26 2019

@author: gvanriessen
"""

import multiprocessing as mp
import time
import os
import imageio
import h5py
import copy
from termcolor import colored

use_gpu = False

if use_gpu:
    import afnumpy as np
    import afnumpy.fft as fft
    use = 'afnumpy/GPU'
else:
    import numpy as np
    import numpy.fft as fft
use = 'numpy/CPU'

#from wpg import srwlpy as srwl
from wpg.srwlib import SRWLOptD, SRWLWfr
from wpg.beamline import Beamline
from wpg.optical_elements import Use_PP
from wpg.wpg_uti_wf import print_mesh
from wpg.useful_code.wfrutils import get_mesh
from wpg.wavefront import Wavefront
from memory_profiler import profile
from extensions.gvrutils import plotWavefront
from extensions.beamlineExt import oeResample
from wpg.optical_elements import Use_PP, Empty


class driftVol():

    def __init__(self, 
                 zrange,
                 propagation_parameters,
                 outputPath, 
                 saveIntensity = True,
                 saveComplex = False,
                 resampleOutput = 1.0, 
                 zoomXOutput = 1.0, 
                 zoomYOutput = 1.0, 
                 method=1):
        """
        param wf is the input wavefield
        param zrange should be a 3-value tuple in order (start Z, stop Z, number of zplanes)
        param propagation_parameters: SRW propagation parameters list
        param crop: crop output at each slice to dimensions in pixels given as tuple (x,y)
        param save:  save nothing (save=0), save h5 full wavefield (svawe=1), save snapshot, optionally cropped (Save=2)
        """     
    
        if not isinstance(zrange, tuple):
            raise TypeError("zrange must be a tuple, (start Z, stop Z, number of zplanes), type [float, float, int]")
        if len(zrange) != 3:
            raise ValueError("zrange tuple must have three values, (start Z, stop Z, number of zplanes), type [float, float, int]")
        if (not isinstance(zrange[0], float)) or (not isinstance(zrange[1], float)) or (not isinstance(zrange[2], int)):
            raise ValueError("zrange tuple must have three values, (start Z, stop Z, number of zplanes), type [float, float, int]")
       
        if outputPath is not None:
            self.outPath = outputPath
            try: 
                os.makedirs(os.path.join(self.outPath,'complex'))
                os.makedirs(os.path.join(self.outPath,'intensity')) 
            except OSError:            
                print('Problem creating output directory')
        else:
            self.outPath = None
            
        self.h5dst = None
            
        self.zrange = zrange
        self.numberZPlanes = zrange[2]
        self.pp = propagation_parameters     
        self.vol = None
        self.data = []   
        self.resampleOutput = resampleOutput
        self.zoomOutput = [zoomXOutput, zoomYOutput]
        self.method=method     
        self.saveComplex   =  saveComplex
        self.saveIntensity = saveIntensity
        self.propagationDistance = 0
        
        print('Progation parameters:' + self.pp.__str__())


    def constructVol(self,wf):
        """
        Try to deduce the output array size based on propagation parameters
        !! Currently assume that autoresize options are not use
        """
        assert self.pp.auto_resize_after ==0 and self.pp.auto_resize_before ==0 , \
               'Not able to predict output wf size when options auto_resize_before or auto_resize_after have value 1'

        nx=int(wf.params.Mesh.nx * self.pp.zoom_h * self.pp.sampling_h)
        ny=int(wf.params.Mesh.ny * self.pp.zoom_v * self.pp.sampling_v)
        nz = self.zrange[2]

        # construct a 3D array with dimensions nx*nx*nz to hold part (real or imag) or output wavefield
        self.vol = np.zeros((nx,ny,nz), dtype=np.float32)

        print('Output dimensions expected: %d x %d x %d' % (nx,ny,nz))

        # nth slice is indexed as vol[:, :, n]


    def propagate(self,wfr, plot=False,outFilePath=None):
    
    # wrapper for compatibility.  
    # OutFilePath changes value of self.outPath
    
        if outFilePath is not None:
            self.outPath =outFilePath
        
        if self.method == 3:
            self.propagateMultiLoky(wfr)
        elif self.method == 2:
            self.propagateMulti(wfr)
        elif self.method == 1:
            self.propagateSequential(wfr)
        else:
            print('Propagation implementation method not recognised.')
        
        # Replaced...
        #now set input wf to the output of final propagation step  
        #wfLastZ= Wavefront()        
        
        #wfLastZ.load_hdf5(os.path.join(os.path.join(self.outPath,'complex'), str(self.zrange[2]-1) + '_' + str(self.zrange[1]) + '.h5')  )
        #wfr = Wavefront(srwl_wavefront=wfLastZ._srwl_wf)
        
  
    def propagateMulti(self,wf):

        starttime = time.time()
        processes = []
        for i,z in zip(range(0,self.zrange[2]), np.linspace(self.zrange[0],self.zrange[1],self.zrange[2])):
            if i == self.zrange[2]-1:
                copyWFR = False
            else:
                copyWFR = True

            p = mp.Process(target=self.propagateToZ, args=(wf,z,i,copyWFR))
            processes.append(p)
         
        [proc.start() for proc in processes]
        
        [proc.join() for proc in processes]

        print('Propagation over {} slices of total \
              length {} completed in {} s.'.format(self.zrange[2], 
                                                   self.propagationDistance, 
                                                   time.time() - starttime))
        
        
        
    def propagateSequential(self,wf):

        starttime = time.time()
        processes = []
        for i,z in zip(range(0,self.zrange[2]), np.linspace(self.zrange[0],self.zrange[1],self.zrange[2])):
            if i == self.zrange[2]-1:
                copyWFR = False
            else:
                copyWFR = True
            
            self.propagateToZ(wf,z,i,copyWFR)
            
        print('Sequential propagation over volume took {} seconds'.format(time.time() - starttime))    
    
                       
    @profile
    def propagateMultiLoky(self,wf):
                
        from loky import get_reusable_executor
        
        # Create an executor with 12 worker processes, that will
        # automatically shutdown after idling for 2s
        executor = get_reusable_executor(max_workers=6, timeout=2)

        starttime = time.time()
        params = []
        indices = range(0,self.zrange[2])
        zVals = np.linspace(self.zrange[0], 
                            self.zrange[1],
                            self.zrange[2])
        
        for i,z in zip(indices, zVals):
            if i == self.zrange[2]-1:
                copyWFR = False
            else:
                copyWFR = True
            params.append( [wf, z, self.outPath, i, copyWFR])
                     
        results = executor.map(self.propagateToZ_params,params)
        
        n_workers = len(set(results))
        
        print(self.data)
    

    def blDrift(self,distance):
    # generate a beamline representing drift region

        bl = Beamline()#description='Drift %s' % str(distance))
        bl.append(SRWLOptD(distance),self.pp)

        return bl
    
    def drift(self,z, wf):        
        d = self.blDrift(z)        
        d.propagate(wf)
        del d
        self.propagationDistance += z


    def resample(self, wf):
        # simplify rescaloing wf
        '''
        if self.zoomOutput != 1 or self.resampleOutput != 1:
            bl = Empty(  )
            bl.propagate(wf, Use_PP(zoom = self.zoomOutput,sampling=self.resampleOutput))        
            del bl
        '''
        bl= Beamline()#description='Drift %s' % str(distance))
        bl.append(SRWLOptD(0),
                  Use_PP(zoom_h = self.zoomOutput[0], 
                         zoom_v = self.zoomOutput[1], 
                         sampling=self.resampleOutput)
                  )
        
        wf_mesh = wf.params.Mesh
        nxi, nyi = wf_mesh.nx, wf_mesh.ny
        
        bl.propagate(wf)
        
        wf_mesh = wf.params.Mesh
        self.nxf, self.nyf = wf_mesh.nx, wf_mesh.ny
        
        print ('Size before resample: {}x{}.  After: {}x{}'.format(nxi,nyi,self.nxf,self.nyf))
        

    def propagateToZ_params(self, params):
    # wrapper for propagateToZ    
       [wf0,z,index,copy] = params
       
       return self.propagateToZ(wf,z,index,copy) 
   
    
    def propagateToZ(self,wf0,z,index,propagateCopy):
       # propagates wavefield over a distance z (not to z!)
       if not isinstance(wf0, Wavefront):
           raise TypeError('First parameter not of type Wavefront')
       
       if propagateCopy==True:
           wf = Wavefront(srwl_wavefront=copy.deepcopy(wf0._srwl_wf))              
       else:
           wf = wf0

       self.drift(z, wf)       
       self.resample(wf)  
       self.save(wf, index, z)
       
       if copy==True: 
           del wf
     
       return index
   
    def save(self, wf, i, z):
       # save the wavefield to file  
       
       baseName =  '{}_{}'.format(i,z)

       if self.outPath is not None and self.saveComplex is True:
           path = os.path.join(os.path.join(self.outPath, 'complex'),baseName+'.h5')
           wf.store_hdf5(path) 
           print ('Wrote wavefield %s' % (path))
           
       if self.outPath is not None and self.saveIntensity is True: 
           
           # block below tested, but need to think about organisation of h5 file for compatibility with 3rd party codes
           if self.h5dst is None:
               self.h5Path = os.path.join(self.outPath, 'intensity3D.h5')
               
               try:
                   self.h5File =  h5py.File(self.h5Path,'w')                 
                   print ('Created {}'.format(self.h5Path))
               except:
                   print ('Failed to create {}'.format(self.h5Path))
               else:
                   self.h5dst =self.h5File.create_dataset("3D Wavefield", 
                                                           shape=(self.numberZPlanes, self.nyf, self.nxf,1),
                                                           dtype=np.uint16)
                   
           if self.h5dst is not None: 
               self.h5dst[i] = wf.get_intensity()
               print ('Wrote wavefield %d to %s' % (i, self.h5Path))
           
           
           path = os.path.join(os.path.join(self.outPath, 'intensity'), f'{i:06}.tif') 
           imageio.imwrite(path, wf.get_intensity() )
           print ('Wrote wavefield %s' % (path))
           

    def showVol(self, resampleFactor=1.0):
        """
           Render 3D wavefield data
           This needs some improvement.  It is recommended that you use TomViz instead
        """
        import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'
        try:
            import visvis as vv
        except ImportError:
            raise ImportError(import_err_msg.format('visvis', 'visvis'))

        if self.vol is None:  #to: correction to use self.data
           
            for i,z in zip(range(0,self.zrange[2]), np.linspace(self.zrange[0],self.zrange[1],self.zrange[2])):

                path = os.path.join( os.path.join(self.outPath, 'complex'), '{}_{}.h5'.format(i,z))

                wfr = Wavefront()
                wfr.load_hdf5(path)
                
                if i==0: 
                    self.constructVol(wfr)
                
                self.vol[:,:,i] = wfr.get_intensity(slice_number=0,polarization='horizontal')

        app = vv.use()

        # Show
        vv.figure()
        #a1 = vv.subplot(121)
        #t1 = vv.volshow(self.vol[:,:,:], renderStyle = 'mip')
        #vv.title('color MIP render')
        #a2 = vv.subplot(122)
        t2 = vv.volshow(self.vol[:,:,:], renderStyle = 'iso')
        t2.isoThreshold = 0.5
        vv.title('color ISO-surface render')

        # Share cameras
        #a1.camera = a2.camera

        # Run app
        app.Run()

def testDriftVol():
    from extensions.blLTZP0 import constructTestWaveField
    from wpg.srwlib import SRWLOptZP, SRWLOptD, SRWLOptA
    from wpg.beamline import Beamline
    from wpg.optical_elements import Use_PP
    #from extensions.beamlineExt import bl

    ekeV = 0.6
    wf=constructTestWaveField(npoints=128,ekeV=ekeV,z1=5)
   
    beamline = Beamline()#(description='bl')
     
    #drift 1 m to ZP
    beamline.append(SRWLOptD(1),Use_PP(semi_analytical_treatment=1.0))
     
    #Central stop defined by an aperture
    beamline.append(SRWLOptA(_shape='c', _ap_or_ob='o', _Dx=10.e-6, _Dy=10.e-6),
                    Use_PP(semi_analytical_treatment=4.0) )

    #Zone plate
    N =400
    dr = 0.1e-06
    beamline.append(SRWLOptZP(_nZones=N, _rn=dr, 
                        _thick=0.5e-06, _delta1=0.0031406, 
                        _atLen1=0.277984e-6), 
               Use_PP())

    wl =  12.4*1e-10/ekeV    
    f =  4*N*(dr**2)/wl
    NA = wl / (2 * dr)
    D  = f*NA
    dof = wl / ( 2 * NA**2 )
    res = 0.610*wl/NA 
    print ('ZP:\n focal length={}\n diameter={}\n NA={}\n dof={}\n resolution={}\n'.format(f, D, NA, dof,res))
    
    zRange = 20*dof
    zPlanes=50
    
    #drift  to  focus - DOF
    beamline.append(SRWLOptD(f-zRange/2),
              Use_PP(semi_analytical_treatment=4.0,zoom=0.3, sampling=10))
    
    #Order sorting aperture 
    beamline.append(SRWLOptA(_shape='c', _ap_or_ob='a', _Dx=9.e-6, _Dy=9.e-6),
                    Use_PP() )

    print(beamline.__str__())
    
    # Propagate through elements defined above
    beamline.propagate(wf)
    
    #plotWavefront(wf, 'Wavefront before drift', cuts=True)
 
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf)
    dx, dy = (xmax-xmin)/nx, (ymax-ymin)/ny
    print('Mesh grid period: {}x{}'.format(dx,dy))
    
    # define driftVol
    distances=(0.0,zRange,zPlanes)
    outPath='/opt/wpg/wpg/extensions/out/driftVol/'
    dVol=driftVol(distances,
                  Use_PP(semi_analytical_treatment=1.0),
                  outPath, 
                  saveIntensity=True, 
                  zoomOutput = 0.1, 
                  method=1)
    
    #beamline.append(dVol,Use_PP())
    
    # Propagate through driftVol
    dVol.propagate(wf)
  
    #dVol.showVol()

def testLTZPFocus():
    
    # test for Alaleh's project
    # propagating input wfr to some drift position and then over some range zrange
    
    from wpg.srwlib import SRWLOptD, SRWLOptA
    from wpg.beamline import Beamline
    from wpg.optical_elements import Use_PP
    from wpg.wavefront import Wavefront
    #from extensions.beamlineExt import bl

    initialWavefront='/opt/wpg/wpg/extensions/in/50Slices63dLTZP_simulationData/wf_noCS.h5'
    outPath='/opt/wpg/wpg/extensions/out/driftVol/kinoform-noCS/'
    
    ekeV = 0.6
    dr = 200.e-9  # correct, Alaleh?    
    wl =  12.4*1e-10/ekeV    
    f =   1.4345e-3 #m   f =  4*N*(dr**2)/wl
    wl =  12.4*1e-10/ekeV     
    NA = wl / (2 * dr)
    D  = f*NA
    dof = wl / ( 2 * NA**2 )
    res = 0.610*wl/NA   #check, guide only
    print (colored('ZP:\n focal length={}\n diameter={}\n NA={}\n dof={}\n resolution={}\n'.format(f, D, NA, dof,res),'red'))
    
    zRange = 300.e-6 #20*dof
    zPlanes= 300
    
    wf = Wavefront()
    wf.load_hdf5(initialWavefront)
    plotWavefront(wf, 'Wavefront input', cuts=True)
    
    beamline = Beamline()#(description='bl')
         
    #Central stop defined by an aperture
    #beamline.append(SRWLOptA(_shape='c', _ap_or_ob='o', _Dx=10.e-6, _Dy=10.e-6),
    #                Use_PP() )
    
    
    #drift  to near focus, specifically f-zRange/2
    beamline.append(SRWLOptD(f-zRange/2),
              Use_PP(semi_analytical_treatment=4.0,sampling=1,zoom_h=0.8, zoom_v=0.5,))  
    #Order sorting aperture 
    #beamline.append(SRWLOptA(_shape='c', _ap_or_ob='a', _Dx=9.e-6, _Dy=9.e-6),
    #                Use_PP() )
    
    #print(beamline.__str__())
    
    # Propagate through elements defined above
    beamline.propagate(wf)
    
    plotWavefront(wf, 'Wavefront before volume generated', cuts=True)
     
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf)
    dx, dy = (xmax-xmin)/nx, (ymax-ymin)/ny
    print('Mesh grid period: {}x{}'.format(dx,dy))
    
    # define driftVol
    distances=(0.0,zRange,zPlanes)
    
    dVol=driftVol(distances,
                  Use_PP(semi_analytical_treatment=1.0,zoom=1, sampling=1),
                  outPath, 
                  zoomXOutput = 0.6,     # reduce size of output, keep all significant features
                  zoomYOutput = 0.14,     # reduce size of output, keep all significant features
                  resampleOutput=0.5,     # reduce size of output, proportional loss of resolution
                  method=1,
                  saveIntensity=True,
                  saveComplex=False)
    
    #beamline.append(dVol,Use_PP())
    
    # Propagate through driftVol
    dVol.propagate(wf)



def testLTZPOverview():
    
    # test for Alaleh's project
    # propagate from input wavefield, which we assume to be the ESW of an optic, 
    # to some distance given by zrange
      
    from wpg.srwlib import SRWLOptD, SRWLOptA
    from wpg.beamline import Beamline
    from wpg.optical_elements import Use_PP
    from wpg.wavefront import Wavefront
    #from extensions.beamlineExt import bl

   # initialWavefront='/opt/wpg/wpg/extensions/in/50Slices63dLTZP_simulationData/wf_noCS.h5'
    initialWavefront='/opt/wpg/wpg/extensions/in/esw_h5Files/wfPostZPX_50Slices_withCS_356nmZ161nmX.h5'
    outPath='/opt/wpg/wpg/extensions/out/driftVol2/50Slices_withCS_356nmZ161nmX/'
    
    ekeV = 0.6
    dr = 200.e-9  # correct, Alaleh?    
    wl =  12.4*1e-10/ekeV    
    f =   1.4345e-3 #m   f =  4*N*(dr**2)/wl
    wl =  12.4*1e-10/ekeV     
    NA = wl / (2 * dr)
    D  = f*NA
    dof = wl / ( 2 * NA**2 )
    res = 0.610*wl/NA   #check, guide only
    print (colored('ZP:\n focal length={}\n diameter={}\n NA={}\n dof={}\n resolution={}\n'.format(f, D, NA, dof,res),'red'))
    
    zRange = f+200.e-6 #20*dof
    zPlanes= 400
    
    
    wf = Wavefront()
    wf.load_hdf5(initialWavefront)
    plotWavefront(wf, 'Wavefront input', cuts=True)
    
    
    
    #beamline = Beamline()#(description='bl')
         
    #Central stop defined by an aperture
    #beamline.append(SRWLOptA(_shape='c', _ap_or_ob='o', _Dx=10.e-6, _Dy=10.e-6),
    #                Use_PP() )
    
    
    #drift  to near focus, specifically f-zRange/2
    #beamline.append(SRWLOptD(f-zRange/2),
    #          Use_PP(semi_analytical_treatment=4.0,sampling=1,zoom_h=0.8, zoom_v=0.5,))  
    #Order sorting aperture 
    #beamline.append(SRWLOptA(_shape='c', _ap_or_ob='a', _Dx=9.e-6, _Dy=9.e-6),
    #                Use_PP() )
    
    #print(beamline.__str__())
    
    # Propagate through elements defined above
    #beamline.propagate(wf)
    
    #plotWavefront(wf, 'Wavefront before volume generated', cuts=True)
     
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf)
    dx, dy = (xmax-xmin)/nx, (ymax-ymin)/ny
    print('Mesh grid period: {}x{}'.format(dx,dy))
    
    # define driftVol
    distances=(0.0,zRange,zPlanes)
    dVol=driftVol(distances,
                  Use_PP(semi_analytical_treatment=1.0,zoom=1, sampling=1),
                  outPath, 
                  zoomXOutput = 0.2,     # reduce size of output, keep all significant features
                  zoomYOutput = 0.25,     # reduce size of output, keep all significant features
                  resampleOutput=0.5,     # reduce size of output, proportional loss of resolution
                  method=1,
                  saveIntensity=True,
                  saveComplex=False)
    
    #beamline.append(dVol,Use_PP())
    
    # Propagate through driftVol
    dVol.propagate(wf)
  

def main():
    #pass
    #testLTZPFocus()
    testLTZPOverview()

if __name__ == '__main__':
    main()

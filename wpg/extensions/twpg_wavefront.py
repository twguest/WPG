# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:17:08 2019

@author: twguest
"""

import array
import warnings
import os
import numpy as np
import h5py
import pylab as plt
import wpg.srwlib as srwlib
from skimage.external import tifffile
import wpg.utils as utils
import wpg.glossary as glossary

from wpg.utils import srw_obj2str
import imageio
from extensions.process import get_rms
from wpg.srwlib import SRWLOptD,SRWLOptA,SRWLOptC,SRWLOptT,SRWLOptL,SRWLOptMirEl

#from extensions.twpg_uti_wf import coherence
#from extensions.twpg_uti_wf import fixPhase

from wpg.wavefront import Wavefront as SRWLWavefront

from extensions.twpg_beamline import Beamline
from extensions.utils.logging import *
from copy import deepcopy

from extensions.xrt_coherence import calc_coherence
class Wavefront(SRWLWavefront):
    
    def __init__(self, srwl_wavefront = None):
        super().__init__(srwl_wavefront)
        
        self.type = 'se'
        self.eField = self.data.arrEhor    
        self.II = None
        
    def rescale(self, scale_x, scale_y):
        """
    
        :param scale_x: pixel magnification along x-axis
        :param scale_y: pixel magnification along y_axis
        """
    
        
        def rescale_param(scale_x, scale_y = None):
            
            if scale_y == None:
                scale_y = scale_x
            ppscale = [0, 0, 1.0, 0, 0, 1, scale_x, 1, scale_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            return ppscale
    
    
            
        bl = Beamline("Temporary Beamline")
        bl.append(SRWLOptD(0), rescale_param(scale_x, scale_y))
        bl.propagate(self)   
        
        del bl
        

    
    def zoom(self, zoom_x, zoom_y = None):
        """
        :param scale_x: magnification along x-axis
        :param scale_y: magnification along y_axis
        """
        if zoom_y == None:
            zoom_y = zoom_x
    
        def zoom_param(zoom_x, zoom_y):
            ppzoom = [0, 0, 1.0, 0, 0, 1/zoom_x, zoom_x, 1/zoom_y, zoom_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            return ppzoom
    
        bl = Beamline("Temporary Beamline")
        bl.append(SRWLOptD(0), zoom_param(zoom_x, zoom_y))
        bl.propagate(self)   
        
        del bl
        
    def save_tif(self, filename):
        imgReal = self.II[:,:,0]
        imageio.imwrite(filename + ".tif", imgReal)
        
        #tifffile.imsave(filename, imgReal, bigtiff = True)
    def electric_field(self):
        
        itr = 0 
        if self.type == 'gsm':
            self.eField = np.zeros([self.modes[0].data.arrEhor.shape[0],
                self.modes[0].data.arrEhor.shape[1],
                self.modes[0].data.arrEhor.shape[2],
                self.modes[0].data.arrEhor.shape[3]])
        else:
            pass

        if self.type == 'se':
            self.eField = self.data.arrEhor
            
        else: 
            for itr in range(len(self.modes)):
                self.eField += self.modes[itr].data.arrEhor*self.eigenval[itr]
            
    
    def intensity(self):

        
        itr = 0
        
        if self.type == 'se':
            

            self.II = 0
            self.II = self.get_intensity()[:,:,:]
            
        elif self.type == 'gsm' or 'me':
            
            self.II = 0
    
            for mode in self.modes:
                self.II += mode.get_intensity()*self.eigenval[itr]
                itr += 1
        
        return self.II
                
    def collapse(self):
        
        """ 
        Collapses each of the gaussian modes into a single wavefront
        """
        
        self.log.debug("Collapsing Wavefront")
        
        self.data = deepcopy(self.modes[0].data)
        self.params = deepcopy(self.modes[0].params)
        
        self.intensity()
        self.electric_field()
         
        fixPhase(self.eField)
        

    def coherence(self, axisName = 'x', plot = True):

        if axisName == 'x':
            axis = np.linspace(self.params.Mesh.xMin, self.params.Mesh.xMax, self.params.Mesh.nx)
        elif axisName == 'y':
            axis = np.linspace(self.params.Mesh.yMin, self.params.Mesh.yMax, self.params.Mesh.ny)

        coh = coherence(self.eField[:,:,:,0], axisName, axis)

        self.log.info("Intensity Width: {}".format(coh[3]**0.5))
        self.log.info("Coherence Length: {}".format(coh[4]**0.5))
        self.log.info("Transverse Degree of Coherence: {}".format(tdoc))

        beta = (coh[4]**0.5)/(coh[3]**0.5)
        self.log.info("Beta: {}".format(beta))
        
        if plot == True:
            fig = plot_coherence(coh, self.II, axisName, axis)
        print(coh) 
        return coh

    
    def save_intens(self, filename):
        ii = self.get_intensity()[:,:,0]
        imageio.imwrite(filename + ".tif", ii)
        
    def plot(self):
       
        try: 
            plt.figure()
            plt.imshow(self.II[:,:,0])
        except(NameError):
            print("Single or Non-Collapsible Wavefront")
        
    def pixelsize(self):
        
        px = (self.params.Mesh.xMax-self.params.Mesh.xMin)/self.params.Mesh.nx
        py = (self.params.Mesh.yMax-self.params.Mesh.yMin)/self.params.Mesh.ny
        
        return px, py
    
    def write(self, filename):
        output = open(filename + ".txt", "w")
        output.write(self.__str__())
        output.close()
        

    def add_wfr(self, wfr2, eigenval):
        self.data.arrEhor += wfr2.data.arrEhor * eigenval
        self.data.arrEver += wfr2.data.arrEver * eigenval
        
        if self.II is None:
            self.II = self.get_intensity()
        else:
            self.II += wfr2.get_intensity() * eigenval
        return self
   
    def add_int(self, wfr2, eigenval):
        self.II = self.get_intensity()
        self.II += wfr2.get_intensity()*eigenval
        
        self.ph = self.get_phase()
        self.ph += wfr2.get_phase()*eigenval
        
        return self
    def center(self):
        px,py = self.pixelsize()
        
        E = self.data.arrEhor[:,:,0,0]

        coords = np.unravel_index(E.argmax(), E.shape)
        xC = coords[0] * px
        yC = coords[1] * py

        xC += self.params.Mesh.xMin
        yC += self.params.Mesh.yMin
        
        print("Center Coordinate: ({} m, {} m)".format(xC, yC))   
        
        return xC, yC
    

        
        

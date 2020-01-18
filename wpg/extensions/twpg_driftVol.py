#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:53:26 2019

@author: twguest
"""

import multiprocessing as mp
import time
import os
import imageio

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
from extensions.twpg_uti_wf import print_mesh
from extensions.twpg_wavefront import Wavefront
from memory_profiler import profile

from extensions.driftVol import driftVol as DV


#
#def build_empty_wavefront_xy(nx, ny):
#    """
#    Build empty wavefront
#
#    :param nx: Number of point along x-axis
#    :param ny: Number of point along y-axis
#    :return: wpg.Wavefront structure
#
#    """
#
#    wfr = SRWLWfr()  # Initial Electric Field Wavefront
#    wfr.allocate(_ne=1, _nx=nx, _ny=ny)
#
#
#    return wfr



class driftVol(DV):

    
    def __init__(self, zrange,propagation_parameters,outputPath, save=2, cropArea=None):
        super().__init__(zrange,propagation_parameters,outputPath, save=2, cropArea=None)
        
    def propagateToZ(self,wf,z,filePath,index):

       drift = self.blDrift(z)
       drift.propagate(wf)

       # save l writing to same HDF file is too difficult to consider for now)
       wf.save_tif(filePath.split(".")[0] + "slice_".format(index))
       print ('Wrote tif file at plane %d to %s' % (index,filePath))

       #from gvrutils import plotWavefront

       # insert the wavefield intensity into a 3D volume... rethink?
       #self.vol[:,:,index] = wf.get_intensity(slice_number=0,polarization='horizontal')
       
       return index, filePath

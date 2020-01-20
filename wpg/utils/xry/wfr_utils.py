# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import numpy
import pylab
import numpy as np    

from extensions.twpg_wavefront import Wavefront

from sklearn.preprocessing import minmax_scale as normalise


def zoom_wfr(wfr, zoom_x, zoom_y = None, show = True):
    """
    Propagate wavefront through beamline.

    :param wfr: Input wavefront (will be re-writed after propagation)
    :type wfr: wpg.wavefront.Wavefront
    """
    if zoom_y == None:
        zoom_y = zoom_x

    if isinstance(wfr, Wavefront):
        wfr = Wavefront(srwl_wavefront=wfr._srwl_wf)
    

        
    bl = Beamline("Temporary Beamline")
    bl.append(SRWLOptD(0), zoom_param(zoom_x, zoom_y))
    bl.propagate(wfr)   
    
    if show == True:
        plot(wfr)
    del bl
    
def rescale_param(scalefactor_x, scalefactor_y = None):
    if scalefactor_y == None:
        scalefactor_y = scalefactor_x
    ppzoom = [0, 0, 1.0, 0, 0, 1, scalefactor_x, 1, scalefactor_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return ppzoom

def rescale_wfr(wfr, rescale_x, rescale_y):
    """
    Propagate wavefront through beamline.

    :param wfr: Input wavefront (will be re-writed after propagation)
    :type wfr: wpg.wavefront.Wavefront
    """

    if isinstance(wfr, Wavefront):
        wfr = Wavefront(srwl_wavefront=wfr._srwl_wf)
    

        
    bl = Beamline("Temporary Beamline")
    bl.append(SRWLOptD(0), rescale_param(rescale_x, rescale_y))
    bl.propagate(wfr)   
    
    del bl

def sum_modes(nmodes, indir, eigenval, outfile, norm = True, save = False):
    
    data = []

    wfr = Wavefront()
    wfr2 = Wavefront()
    for itr in range(0,nmodes):
        
        if norm == True:
            if itr == 0:
                wfr.load_hdf5(indir + "wfr_mode_{}.hdf5".format(itr))
                #print(wfr.data.arrEhor[:,:,0,0])
                wfr.data.arrEhor[:,:,0,0] = normalise(wfr.data.arrEhor[:,:,0,0])
            else:
                print("Adding Mode: {}".format(itr))
                wfr2.load_hdf5(indir + "wfr_mode_{}.hdf5".format(itr))
                wfr2.data.arrEhor[:,:,0,0] = normalise(wfr2.data.arrEhor[:,:,0,0])
                wfr.add_wfr(wfr2,eigenval[itr])
        
        else:
            if itr == 0:
                wfr.load_hdf5(indir + "wfr_mode_{}.hdf5".format(itr))
                wfr.data.arrEhor = wfr.data.arrEhor
            else:
                print("Adding Mode: {}".format(itr))
                wfr2.load_hdf5(indir + "wfr_mode_{}.hdf5".format(itr))            
                wfr = wfr.add_wfr(wfr2,eigenval[itr])

        if save == True:
            wfr.store_hdf5(outfile + "wfr_product_{}".format(itr))
        
        [Ix, Jx], [Wx, Wy], beta = wfr.get_coherence()
        data.append([Wx, Wy, beta])     
    data = np.asarray(data)
    np.savetxt(outfile + "coherence.csv",data,delimiter = ",")
    return wfr

def ff_sum(wfr1, wfr2):
    wfr1.II = wfr1.get_intensity()
    wfr2.II = wfr2.get_intensity()
    
    wfr1.ph = wfr1.get_phase()
    wfr2.ph = wfr2.get_phase()
    ffwfr = Wavefront()
    ffwfr.II = np.zeros(wfr1.II.shape, wfr1.II.dtype)
    
    ffwfr.II += wfr1.II + wfr2.II + 2*np.sqrt(wfr1.II*wfr2.II)*np.cos(wfr1.ph-wfr2.ph)
    return ffwfr


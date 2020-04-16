#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:51:10 2020

@author: twguest
"""

import scipy

import numpy as np

from copy import deepcopy

from wpg import srwlpy

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import calculate_fwhm

from scipy.constants import c

h = scipy.constants.physical_constants['Planck constant in eV s'][0]

def get_wavelength(wfr):
    return (h*c)/(wfr.params.photonEnergy)


def toKspace(wfr):
    """
    Changes the representation of wfr from space-time to space-freq
    
    :param wfr: wpg wavefront strucure in space-time domain
    :returns kfr: wpg wavefront structure in space-freq domain
    """    
    srwl_a = deepcopy(wfr._srwl_wf)
    srwl_a = srwlpy.SetRepresElecField(srwl_a, 'a')
    kfr = Wavefront(srwl_a)
    
    return kfr

def calcDivergence(wfr):
    """
    calculate the full-angle divergence of the beam
    
    :param wfr: WPG wavefront structure
    """
    
    kfr = toKspace(wfr)
    
    return [calculate_fwhm(kfr)['fwhm_x'], calculate_fwhm(kfr)['fwhm_y']]

def sampling(wavelength, propD, roi, dx_i):
    """
    :param wavelength: source wavelength [m]
    :param propD: propagation distance [m]
    :param roi: region of interest (extent of wavefield)
    :param pixel_i: current pixel size [m]
    :returns pixel_f: pixel size in propagated plane [m]
    """
    dx_f = (wavelength*propD)/(roi*dx_i*np.sqrt(2))
    
    return dx_f
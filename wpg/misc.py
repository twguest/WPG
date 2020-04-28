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

def fresnel_sampling(z, wav, dx1, D2, D1):
    """
    :param z: propagation distance
    :param wav: source wavelength
    :param dx1: pixel size in unpropagated plane
    :param D2: width of propagated plane
    :param D1: width of unpropagated plane
    """
    dx2 = (wav*z+dx1*D2)/(D1)
    return dx2
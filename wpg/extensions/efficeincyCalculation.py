#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 14:37:09 2019

@author: Al
"""
#Calculate efficiency 

from wpg.optical_elements import _save_object, _load_object, Empty
from wpg.srwlib import srwl_uti_ph_en_conv
from wpg.generators import build_gauss_wavefront_xy
from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, plot_intensity_map, averaged_intensity,calculate_fwhm, get_intensity_on_axis
from wpg.useful_code.wfrutils import propagate_wavefront, calculate_fwhm_y,calculate_fwhm_x, print_beamline, get_mesh

import os, sys, imageio
from slugify import slugify

from extensions.gvrutils import propagationParameters, calculate_theta_fwhm_cdr, \
                                plotWavefront,plot_wfront2, blSave, bl
from extensions.gvrutils import writeIntensity, pixelScale
from extensions.gvrutils import get_transmissionFunctionCplx, calculate_fwhm
from extensions.gvrutils import get_intensity_on_axis, printWFProperties
from extensions.gvrutils import show_transmissionFunctionCplx,  visVol,Struct, blDrift, blRescale, BeamlineCustom
from extensions.multisliceOptE import *


# Read h5 file at the focus
atFocus = '/opt/wpg/wpg/extensions/out/2019-05-23/400eb000-7d2f-11e9-8965-194a49d9ce2f/wfFocusX.h5'

# Read h5 file after CS, immidiately before LTZP


# Calculate total I at First Order diffraction
wf = Wavefront()
wf.load_hdf5(atFocus)

firstOrderI = wf.get_intensity().sum()
print ('First Order Intensity =', firstOrderI)
print ('FWHM [um]:', calculate_fwhm_y(wf) *1.e6)

writeIntensity(wf,
              section.description,
              path=strOutPath,
              polarization='horizontal',
              imgType='tif')
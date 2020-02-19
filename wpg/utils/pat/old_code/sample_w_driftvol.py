#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:28:32 2019

@author: patrick
"""


import os, sys

sys.path.append(r"WPG")
sys.path.append(r"/opt/wpg/wpg")
sys.path.append(r"/opt/wpg/wpg/extensions")

import matplotlib.pyplot as plt
import imageio
import numpy as np

import shelve
from wpg import (
    srwl_bl,
    srwlib,
    srwl_uti_smp,
    uti_io,
    srwl_uti_src,
    wavefront,
    optical_elements,
    beamline,
    generators,
    wpg_uti_wf,
    srwlpy,
)


linux = True


# open shelf for saved varibles and open them
if linux:
    d = shelve.open("shelve_data/asp_source_w_kb_optics_v4")
else:
    d = shelve.open("shelve_data\\asp_source_w_kb_optics_v4")
op = d["op"]
wf = d["sample_wf"]

# Make a wpg type wavefront
wf_wpg = wavefront.Wavefront(wf)

# save this initial intensity
wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity
wf_wpg_phase = wf_wpg.get_phase(slice_number=0, polarization="total")


if linux:
    imageio.imwrite("io/sample/norm/init_wf.tiff", wf_wpg_inten)
    imageio.imwrite("io/sample/log/log10_init_wf.tiff", np.log10(wf_wpg_inten))
    imageio.imwrite("io/sample/phase/phase_smp_plane.tiff", wf_wpg_phase)
else:
    imageio.imwrite("io\\sample\\norm\\init_wf.tiff", wf_wpg_inten)
    imageio.imwrite("io\\sample\\log\\log10_init_wf.tiff", np.log10(wf_wpg_inten))


from driftVol import driftVol

distances = (0.0, 3.5, 6)
pp = [0, 0, 1.0, 0, 3, 1.3, 0.76, 1.3, 0.76, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
outpath = "/tmp/"

dvol = driftVol(distances, pp, "")

dvol.propagate(wf_wpg, outFilePath="")

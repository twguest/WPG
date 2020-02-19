import os, sys

sys.path.append(r"C:\Users\patri\Documents\Uni\thesis\python\WPG_simulations\WPG")
import matplotlib.pyplot as plt
import numpy as np
import time as t
import scipy.ndimage


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

import array as ar


d = shelve.open("shelve_data\\asp_source_w_kb_optics_v4")
op = d["op"]

wf2 = d["sample_wf"]

wf2_wpg = wavefront.Wavefront(wf2)


px_size = (abs(wf2_wpg._srwl_wf.mesh.xStart) + abs(wf2_wpg._srwl_wf.mesh.xFin)) / abs(
    wf2_wpg._srwl_wf.mesh.nx
)

sample1 = srwl_uti_smp.srwl_opt_setup_transm_from_file(
    "cube_samples.tiff", px_size / 1000, 0.1e-6, 1, 0.1e-6
)
# drift1 = srwlib.SRWLOptD(_L=0.0)

drift2 = srwlib.SRWLOptD(_L=0.0)


# resize_bl = beamline.Beamline()
# resize_bl.append(drift1, [0, 0, 1.0, 0, 0,         .5, 2, .5, 2,         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
# resize_bl.propagate(wf2_wpg)

wf2_wpg_inten = wf2_wpg.get_intensity(slice_number=0)
plt.figure()
plt.imshow(wf2_wpg_inten)


sample_bl = beamline.Beamline()

sample_bl.append(
    sample1, [0, 0, 1.0, 0, 0, 1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)
sample_bl.propagate(wf2_wpg)


detector_bl = beamline.Beamline()
detector_bl.append(
    drift2, [0, 0, 1.0, 0, 0, 5, 0.2, 5, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)
detector_bl.propagate(wf2_wpg)

wf2_wpg_inten = wf2_wpg.get_intensity(slice_number=0)
plt.figure()
plt.imshow(wf2_wpg_inten)

plt.show()

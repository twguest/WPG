
# sample.py
# shift and center the beam after propagation through the beamline,
# import a tiff file transmission function and propagte it through the sample.

import os, sys
sys.path.append(r'WPG')
sys.path.append(r'/opt/wpg/wpg')

import matplotlib.pyplot as plt
import imageio
import numpy as np

import shelve
from wpg import srwl_bl, \
    srwlib, \
    srwl_uti_smp, \
    uti_io, \
    srwl_uti_src, \
    wavefront,\
    optical_elements, \
    beamline, \
    generators, \
    wpg_uti_wf, \
    srwlpy


linux = False




# open shelf for saved varibles and open them
if linux:
    d = shelve.open('shelve_data/asp_source_w_kb_optics_v4')
else:
    d = shelve.open('shelve_data\\asp_source_w_kb_optics_v4')
op = d['op']
wf = d['sample_wf']

# Make a wpg type wavefront
wf_wpg = wavefront.Wavefront(wf)

#save this initial intensity
wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

if linux:
    imageio.imwrite('io/sample/norm/init_wf.tiff', wf_wpg_inten)
    imageio.imwrite('io/sample/log/log10_init_wf.tiff', np.log10(wf_wpg_inten))
else:
    imageio.imwrite('io\\sample\\norm\\init_wf.tiff', wf_wpg_inten)
    imageio.imwrite('io\\sample\\log\\log10_init_wf.tiff', np.log10(wf_wpg_inten))


#get the pixel size of the wavefront
px_size = (abs(wf_wpg._srwl_wf.mesh.xStart) + abs(wf_wpg._srwl_wf.mesh.xFin))/(abs(wf_wpg._srwl_wf.mesh.nx))

#make a transmission object from the tiff
##Magnetite Fe3O4, E= 8.26keV density = 5.17g/cm3, delta=?, att = 9.4e-6m
sample1 = srwl_uti_smp.srwl_opt_setup_transm_from_file('slice_chunk_0_10.tiff', px_size/2, px_size*1000, 1, 9.4e-6, xc=-0e-6,yc=0e-6,extTr = 1)
                                                        ##file_path, resolution, thickness, delta, atten_len, center x,y, outside_val


sample_bl = beamline.Beamline()
sample_bl.append(sample1, [0, 0, 1.0, 0, 0,    1, 1, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
sample_bl.propagate(wf_wpg)



wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

if linux:
    imageio.imwrite('io/sample/norm/smp_plane.tiff', wf_wpg_inten)
    imageio.imwrite('io/sample/log/log10_smp_plane.tiff', np.log10(wf_wpg_inten))
else:
    imageio.imwrite('io\\sample\\norm\\smp_plane.tiff', wf_wpg_inten)
    imageio.imwrite('io\\sample\\log\\log10_smp_plane.tiff', np.log10(wf_wpg_inten))




#loop through propagation distances

for i in range(50):

    print('\n\n')
    print(i, '/50')
    detector_bl = beamline. Beamline()
    drift1 = srwlib.SRWLOptD(_L=0.0001)
    detector_bl.append(drift1, [0, 0, 1.0, 0, 0, 1, 1,1,1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #detector_bl.append(drift1, [0, 0, 1.0, 0, 0,    1.5, 0.666, 1.5, .666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    detector_bl.propagate(wf_wpg)
    wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

    if linux:
        imageio.imwrite('io/sample/norm/detector' + str(i + 1) + '.tiff', wf_wpg_inten)
        imageio.imwrite('io/sample/log/log10_detector' + str(i + 1) + '.tiff', np.log10(wf_wpg_inten))
    else:
        imageio.imwrite('io\\sample\\norm\\detector'+str(i+1)+'.tiff', wf_wpg_inten))
        imageio.imwrite('io\\sample\\log\\log10_detector'+str(i+1)+'.tiff', np.log10(wf_wpg_inten))




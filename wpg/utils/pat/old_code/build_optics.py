import os, sys

sys.path.append("/opt/wpg/wpg")
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import time as t
import shelve
from scipy.misc import imsave
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

start_time = t.time()  # Start a timer

d = shelve.open("shelve_data/asp_source_w_kb_optics_v4")  # load variables from shelf
op = d["op"]  # load optical elements defined in asp_source....
wfs = d["wfs"]  # load wavefronts from build source
pp = d["pp"]  # Load propagation params from build_pp

wf1 = wfs[5]  # select a single wavefront (check 'pixels' list in build_source)

bls = []  # Init list of beamlines, one for each optical element
total_beamline = (
    beamline.Beamline()
)  # Init an empty beamline for all the optical elements

for i in range(
    len(op.arOpt)
):  # loop through each optical element (-2 to chop outlast two elements)
    total_beamline.append(op.arOpt[i], pp[i])  # add the element to the total beamline
    bl = beamline.Beamline()  # make a seperate beamline
    bl.append(op.arOpt[i], pp[i])  # append just the single element to the new beamline
    bls.append(bl)  # append the beamline with a single element to the list of beamlines

wf_wpg = wavefront.Wavefront(wf1)  # Make a WPG type wavefront

# init a list of values to track the array size (detector length) and number of pixels
arr_sizes_x = []
arr_sizes_y = []
num_pxs_x = []
num_pxs_y = []


print("\n Initial Wavefront")
# get array data and append to a list
arr_sizes_x.append(abs(wf_wpg._srwl_wf.mesh.xStart - wf_wpg._srwl_wf.mesh.xFin) * 1)
arr_sizes_y.append(abs(wf_wpg._srwl_wf.mesh.yStart - wf_wpg._srwl_wf.mesh.yFin) * 1)
num_pxs_x.append(wf_wpg._srwl_wf.mesh.nx)
num_pxs_y.append(wf_wpg._srwl_wf.mesh.ny)

wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

# save the intensity and log of intensity
print("Saving Image...")

imsave("io/" + str(0) + ".tiff", wf_wpg_inten)
imsave("io/" + "log10_" + str(0) + ".tiff", np.log10(wf_wpg_inten))


# Propagate the wave field through every optical element
for i in range(len(bls)):
    print("\nOptical element:" + str(i + 1))
    bls[i].propagate(wf_wpg)  # propagate though the current optical element.

    # get array data and append to a list
    arr_sizes_x.append(abs(wf_wpg._srwl_wf.mesh.xStart - wf_wpg._srwl_wf.mesh.xFin) * 1)
    arr_sizes_y.append(abs(wf_wpg._srwl_wf.mesh.yStart - wf_wpg._srwl_wf.mesh.yFin) * 1)
    num_pxs_x.append(wf_wpg._srwl_wf.mesh.nx)
    num_pxs_y.append(wf_wpg._srwl_wf.mesh.ny)

    wf_wpg_inten = wf_wpg.get_intensity(slice_number=0)  # get intensity

    # save the intensity and log of intensity
    print("Saving Image...")
    # imageio.imwrite('io\\' + str(i+1) + '.tiff', wf_wpg_inten)
    # log10_inten = np.log10(wf_wpg_inten)
    # imageio.imwrite('io\\' + 'log10_' + str(i+1) + '.tiff', log10_inten)

    imsave("io/" + str(i + 1) + ".tiff", wf_wpg_inten)
    log10_inten = np.log10(wf_wpg_inten)
    imsave("io/" + "log10_" + str(i + 1) + ".tiff", log10_inten)


d["bls"] = bls
d["total_beamline"] = total_beamline

d["sample_wf"] = wf_wpg._srwl_wf

end_time = t.time()

print("Time taken: " + str(end_time - start_time) + " seconds\n\n")

d.close()


plt.figure()
plt.imshow(wf_wpg_inten)


plt.figure()
plt.plot(arr_sizes_x, ".")
plt.plot(arr_sizes_y, "r.")

plt.figure()
plt.plot(num_pxs_x, ".")
plt.plot(num_pxs_y, "r.")
plt.show()

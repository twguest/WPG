# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os,shutil
import copy
import numpy
import pylab

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

from wpg.srwlib import srwl
from wpg.wavefront import Wavefront

__author__ = "A. Buzmakov, L. Samoylova, C. Fortmann-Grote"


def print_mesh(wf):
    """
    #print out wfr wavefront mesh.

    :param wf: wpg.Wavefront structure
    """

    wf_mesh = wf.params.Mesh
    w_space = wf.params.wSpace
    ##print(w_space)
    if w_space == "R-space":
        print(
            "nx {:5d}  range_x [{:.1e}, {:.1e}] mm".format(
                wf_mesh.nx, wf_mesh.xMin * 1e3, wf_mesh.xMax * 1e3
            )
        )
        print(
            "ny {:5d}  range_y [{:.1e}, {:.1e}] mm".format(
                wf_mesh.ny, wf_mesh.yMin * 1e3, wf_mesh.yMax * 1e3
            )
        )
    if w_space == "Q-space":
        print(
            "nx {:5d}  range_x [{:.1e}, {:.1e}] mrad".format(
                wf_mesh.nx, wf_mesh.qxMin * 1e3, wf_mesh.qxMax * 1e3
            )
        )
        print(
            "ny {:5d}  range_y [{:.1e}, {:.1e}] mrad".format(
                wf_mesh.ny, wf_mesh.qyMin * 1e3, wf_mesh.qyMax * 1e3
            )
        )
    return


def calc_pulse_energy(wf):
    """
    calculate energy of  in time domain

    :param wf: wpg.Wavefront structure
    :return: pulse energy value, J
    """
    J2eV = 6.24150934e18
    if wf.params.wDomain != "time":
        print(
            "Pulse energy cannot be calculated for {:s} domain".format(
                wf.params.wDomain
            )
        )
        return None
    else:
        dx = (wf.params.Mesh.xMax - wf.params.Mesh.xMin) / (wf.params.Mesh.nx - 1)
        dy = (wf.params.Mesh.yMax - wf.params.Mesh.yMin) / (wf.params.Mesh.ny - 1)
        dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / (
            wf.params.Mesh.nSlices - 1
        )
        pulse_energy = wf.get_intensity().sum(axis=0).sum(axis=0).sum(axis=0)
        pulse_energy_J = pulse_energy * dx * dy * 1e6 * dt
        print(
            "Number of photons per pulse: {:e}".format(
                pulse_energy_J * J2eV / wf.params.photonEnergy
            )
        )
        return pulse_energy_J


def averaged_intensity(wf, bPlot=False):
    """
    wrapper for integral_intensity() for backward compatibility

    """
    integral_intensity(wf, bPlot)


def integral_intensity(wf, threshold=0.01, bPlot=True):
    """
    plot the slice-to-slice integral intensity averaged over a meaningful range

    :param wf: wavefront structure
    :param threshold: defined the threshold for slices, integrated_slice_intensity_max*threshold
    :param bPlot: if True plot temporary structure or spectrum in the meaningful range
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold, mainly needed for processing spiky FEL source

    """
    J2eV = 6.24150934e18
    # total0=wf.get_intensity().sum();
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)  # I(slice_num)
    int0 = int0 * (dx * dy * 1.0e6)  # wf amplitude units sqrt(W/mm^2)

    # Get center pixel numbers.
    center_nx = int(mesh.nx / 2)
    center_ny = int(mesh.ny / 2)
    int0_00 = wf.get_intensity()[center_ny, center_nx, :]
    int0max = max(int0)

    # Get meaningful slices.
    aw = [a[0] for a in numpy.argwhere(int0 > int0max * threshold)]
    int0_mean = int0[min(aw) : max(aw)]  # meaningful range of pulse
    if bPlot:
        if mesh.nSlices > 1:
            dSlice = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        else:
            dSlice = 0
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices) * dSlice + mesh.sliceMin, int0)
        pylab.plot(
            numpy.arange(min(aw), max(aw)) * dSlice + mesh.sliceMin, int0_mean, "ro"
        )
        if wf.params.wDomain == "time":
            pylab.title("Power")
            pylab.xlabel("s")
            pylab.ylabel("W")
        else:  # frequency domain
            pylab.title("Spectral Energy")
            pylab.xlabel("eV")
            pylab.ylabel("J/eV")
        pylab.show()
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices) * dSlice + mesh.sliceMin, int0_00)
        pylab.plot(
            numpy.arange(min(aw), max(aw)) * dSlice + mesh.sliceMin,
            int0_00[min(aw) : max(aw)],
            "ro",
        )
        if wf.params.wDomain == "time":
            pylab.title("On-Axis Power Density")
            pylab.xlabel("s")
            pylab.ylabel("W/mm^2")
        else:  # frequency domain
            pylab.title("On-Axis Spectral Fluence")
            pylab.xlabel("eV")
            pylab.ylabel("J/eV/mm^2")
        pylab.show()
    averaged = int0_mean.sum() / len(int0_mean)
    #print("number of meaningful slices:", len(int0_mean))
    if wf.params.wDomain == "time":
        dt = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        #print("Pulse energy {:1.2g} J".format(int0_mean.sum() * dt))
    return averaged


def plot_t_wf(wf, save="", range_x=None, range_y=None, im_aspect="equal"):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map()
    for backward compatibility

    """
    integral_intensity(wf, bPlot=True)
    plot_intensity_map(wf, save, range_x, range_y, im_aspect)


def plot_wf(wf, save="", range_x=None, range_y=None, im_aspect="equal"):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map()
    for backward compatibility

    """
    integral_intensity(wf, bPlot=True)
    plot_intensity_map(wf, save, range_x, range_y, im_aspect)


def plot_intensity_map(wf, save="", range_x=None, range_y=None, im_aspect="equal"):
    """
    Plot wavefront in  R-space.

    :param wf: wavefront structure
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
    :param im_aspect: aspect for 2D image, string or float number, see matplotlib set_aspect().
    """
    import matplotlib.pyplot as plt

    # Get the wavefront and integrate over time.
    wf_intensity = wf.get_intensity().sum(axis=-1)

    # Get average and time slicing.
    # average = averaged_intensity(wf, bPlot=True)
    nslices = wf.params.Mesh.nSlices
    if nslices > 1:
        dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / (nslices - 1)
        t0 = dt * nslices / 2 + wf.params.Mesh.sliceMin
    else:
        t0 = (wf.params.Mesh.sliceMax + wf.params.Mesh.sliceMin) / 2

    # Setup a figure.
    plt.figure(figsize=(10, 10), dpi=100)
    plt.axis("tight")
    # Profile plot.
    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    # Get limits.
    xmin, xmax, ymax, ymin = wf.get_limits()

    # Plot profile as 2D colorcoded map.
    profile.imshow(
        wf_intensity,
        extent=[xmin * 1.0e3, xmax * 1.0e3, ymax * 1.0e3, ymin * 1.0e3],
        cmap="YlGnBu_r",
    )
    profile.set_aspect(im_aspect)

    # Get x and y ranges.
    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    # x = numpy.linspace(xmin*1.e3,xmax*1.e3,wf_intensity.shape[0])
    # y = numpy.linspace(ymin*1.e3,ymax*1.e3,wf_intensity.shape[1])
    x = numpy.linspace(xmin * 1.0e3, xmax * 1.0e3, wf_intensity.shape[1])
    y = numpy.linspace(ymin * 1.0e3, ymax * 1.0e3, wf_intensity.shape[0])

    # Labels.
    profile.set_xlabel("$mm$", fontsize=12)
    profile.set_ylabel("$mm$", fontsize=12)

    # x-projection plots above main plot.
    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    #print(x.shape, wf_intensity.sum(axis=0).shape)
    fwhm = calculate_fwhm(wf)

    #print("FWHM in x = %4.3e m." % (fwhm["fwhm_x"]))
    #print("FWHM in y = %4.3e m." % (fwhm["fwhm_y"]))

    x_projection.plot(x, wf_intensity.sum(axis=0), label="x projection")

    # Set range according to input.
    if range_x is None:
        profile.set_xlim([xmin * 1.0e3, xmax * 1.0e3])
    else:
        profile.set_xlim([-range_x / 2.0, range_x / 2.0])

    # Set title.
    if wf.params.wDomain == "time":
        x_projection.set_title("t0={:03.1g} s ".format(t0))
    else:  # frequency domain
        x_projection.set_title("E0={:05.2g} eV".format(t0))

    # y-projection plot right of main plot.
    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label="y projection")

    # Hide minor tick labels, they disturb here.
    plt.minorticks_off()

    # Set range according to input.
    if range_y is None:
        profile.set_ylim([ymin * 1.0e3, ymax * 1.0e3])
    else:
        profile.set_ylim([-range_y / 2.0, range_y / 2.0])

    # If requested, save to disk, otherwise show in interactive window.
    if save != "":
        # Add parameters.
        plt.savefig(save)
    else:
        plt.show()


# def plot_t_wf_a(wf, save='', range_x=None, range_y=None):
#     """
#     Plot wavefront in Q-space.

#     :param wf: wavefront structure
#     :params: save: Whether to save the figure on disk
#     :type:  string for filename. Empty string '' means don't save.
#     :default: '', do not save the figure.

#     :params: range_x: x-axis range.
#     :type: float
#     :default: None, take entire x range.

#     :params: range_y: y-ayis range.
#     :type: float
#     :default: None, take entire y range.


#     """


def look_at_q_space(wf, output_file=None, save="", range_x=None, range_y=None):
    """
    a wrapper for backward compatibility

    """
    plot_intensity_qmap(wf, output_file, save, range_x, range_y)


def plot_intensity_qmap(
    wf, output_file=None, save="", range_x=None, range_y=None, im_aspect="equal"
):
    """
    change wavefront representation from R- to Q-space and plot it the resulting wavefront.

    :param wf: Wavefront object in R-space representation
    :param output_file: if parameter present - store wavefront in Q-space to the file
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
    :return: propagated wavefront object:
    """
    wfr = Wavefront(srwl_wavefront=wf._srwl_wf)

    if not wf.params.wSpace == "R-space":
        #print("space should be in R-space, but not " + wf.params.wSpace)
        return
    srwl_wf = wfr._srwl_wf
    srwl_wf_a = copy.deepcopy(srwl_wf)
    srwl.SetRepresElecField(srwl_wf_a, "a")
    wf_a = Wavefront(srwl_wf_a)
    if output_file is not None:
        #print("store wavefront to HDF5 file: " + output_file + "...")
        wf_a.store_hdf5(output_file)
        #print("done")

    #print(calculate_fwhm(wf_a))

    # plot_t_wf_a(wf_a, save=save, range_x=range_x, range_y=range_y)
    import matplotlib.pyplot as plt

    wf_intensity = wf_a.get_intensity().sum(axis=-1)
    plt.figure(figsize=(10, 10), dpi=100)
    plt.axis("tight")

    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    xmin, xmax, ymax, ymin = wf_a.get_limits()

    profile.imshow(
        wf_intensity,
        extent=[xmin * 1.0e6, xmax * 1.0e6, ymax * 1.0e6, ymin * 1.0e6],
        cmap="YlGnBu_r",
    )
    profile.set_aspect(im_aspect)

    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    # x = numpy.linspace(xmin*1.e6,xmax*1.e6,wf_intensity.shape[0])
    # y = numpy.linspace(ymin*1.e6,ymax*1.e6,wf_intensity.shape[1])
    x = numpy.linspace(xmin * 1.0e6, xmax * 1.0e6, wf_intensity.shape[1])
    y = numpy.linspace(ymin * 1.0e6, ymax * 1.0e6, wf_intensity.shape[0])
    profile.set_xlabel(r"$\mu$rad", fontsize=12)
    profile.set_ylabel(r"$\mu$rad", fontsize=12)

    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    #print(x.shape, wf_intensity.sum(axis=0).shape)
    x_projection.plot(x, wf_intensity.sum(axis=0), label="x projection")
    if range_x is None:
        profile.set_xlim([xmin * 1.0e6, xmax * 1.0e6])
    else:
        profile.set_xlim([-range_x / 2.0, range_x / 2.0])

    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label="y projection")

    # Hide minor tick labels.
    plt.minorticks_off()

    if range_y is None:
        profile.set_ylim([ymin * 1.0e6, ymax * 1.0e6])
    else:
        profile.set_ylim([-range_y / 2.0, range_y / 2.0])

    if save != "":
        plt.savefig(save)
    else:
        plt.show()

    return


def propagate_wavefront(wavefront, beamline, output_file=None):
    """
    Propagate wavefront and store it in output file.

    :param wavefront: wpg.Wavefront object or path to HDF5 file
    :param beamline: SRWLOptC container of beamline
    :param output_file: if parameter present - store propagaed wavefront to file
    :return: propagated wavefront object
    """

    if not isinstance(beamline, Beamline):
        bl = Beamline(beamline)
    else:
        bl = beamline
    #print(bl)

    if isinstance(wavefront, Wavefront):
        wfr = Wavefront(srwl_wavefront=wavefront._srwl_wf)
    else:
        #print("*****reading wavefront from h5 file...")
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)

    #print_mesh(wfr)
    #print("*****propagating wavefront (with resizing)...")
    bl.propagate(wfr)

    if output_file is not None:
        #print("save hdf5:", output_file)
        wfr.store_hdf5(output_file)
    #print("done")
    return wfr


def calculate_fwhm(wfr):
    """
    Calculate FWHM of the beam calculating number of point bigger then max / 2 throuhgt center of the image

    :param wfr:  wavefront
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m]
    """
    #    intens = wfr.get_intensity(polarization='total')
    intens = wfr.get_intensity(polarization="total").sum(axis=-1)

    mesh = wfr.params.Mesh
    if wfr.params.wSpace == "R-space":
        dx = (mesh.xMax - mesh.xMin) / mesh.nx
        dy = (mesh.yMax - mesh.yMin) / mesh.ny
    elif wfr.params.wSpace == "Q-space":
        dx = (mesh.qxMax - mesh.qxMin) / mesh.nx
        dy = (mesh.qyMax - mesh.qyMin) / mesh.ny
    else:
        return

    x_center = intens[intens.shape[0] // 2, :]
    fwhm_x = len(x_center[x_center > x_center.max() / 2]) * dx

    y_center = intens[:, intens.shape[1] // 2]
    fwhm_y = len(y_center[y_center > y_center.max() / 2]) * dy
    if wfr.params.wSpace == "Q-space":
        #print(wfr.params.wSpace)
        wl = 12.398 * 1e-10 / (wfr.params.photonEnergy * 1e-3)  # WaveLength
        fwhm_x = fwhm_x * wl
        fwhm_y = fwhm_y * wl

    return {"fwhm_x": fwhm_x, "fwhm_y": fwhm_y}


def get_intensity_on_axis(wfr):
    """
    Calculate intensity (e.g. spectrum in frequency domain) along (x=y=0)

    :param wfr:  wavefront
    :return: [z,s0] in [a.u.]
    """

    wf_intensity = wfr.get_intensity(polarization="horizontal")
    mesh = wfr.params.Mesh
    # array dimensions # <-to avoid wrong dimension assignment
    dim = numpy.shape(wf_intensity)
    sz = numpy.zeros(shape=(mesh.nSlices, 2), dtype="float64")
    sz[:, 0] = numpy.linspace(mesh.sliceMin, mesh.sliceMax, mesh.nSlices)
    # <-to avoid wrong dimension assignment
    sz[:, 1] = wf_intensity[dim[0] // 2, dim[1] // 2, :] / wf_intensity.max()

    return sz


def check_sampling(wavefront):
    """ Utility to check the wavefront sampling. """
    xMin = wavefront.params.Mesh.xMin
    xMax = wavefront.params.Mesh.xMax
    nx = wavefront.params.Mesh.nx
    yMin = wavefront.params.Mesh.yMin
    yMax = wavefront.params.Mesh.yMax
    ny = wavefront.params.Mesh.ny
    dx = (xMax - xMin) / (nx - 1)
    dy = (yMax - yMin) / (ny - 1)
    xx = calculate_fwhm(wavefront)
    fwhm_x = xx["fwhm_x"]
    fwhm_y = xx["fwhm_y"]
    Rx = wavefront.params.Rx
    Ry = wavefront.params.Ry
    ekev = wavefront.params.photonEnergy * 1e-3
    dr_ext_x = 12.39e-10 / ekev * Rx / (2 * fwhm_x)
    dr_ext_y = 12.39e-10 / ekev * Ry / (2 * fwhm_y)

    format_string = "|{:4.3e}|{:4.3e}|{:4.3e}|{:4.3e}|{:4.3e}|{:4.3e}|{:4.3e}|"

    ret = "WAVEFRONT SAMPLING REPORT\n"
    ret += "+----------+---------+---------+---------+---------+---------+---------+---------+\n"
    ret += "|x/y       |FWHM     |px       |ROI      |R        |Fzone    |px*7     |px*10    |\n"
    ret += "+----------+---------+---------+---------+---------+---------+---------+---------+\n"
    ret += (
        "|Horizontal"
        + format_string.format(fwhm_x, dx, (xMax - xMin), Rx, dr_ext_x, dx * 7, dx * 10)
        + "\n"
    )
    ret += (
        "|Vertical  "
        + format_string.format(fwhm_y, dy, (yMax - yMin), Ry, dr_ext_y, dy * 7, dy * 10)
        + "\n"
    )
    ret += "+----------+---------+---------+---------+---------+---------+---------+---------+\n\n"

    if 7 * dx < dr_ext_x and 10 * dx > dr_ext_x:
        ret += "Horizontal Fresnel zone extension within [7,10]*pixel_width -> OK\n"
    else:
        ret += 'Horizontal Fresnel zone extension NOT within [7,10]*pixel_width -> Check pixel width."\n'

    if 7 * dy < dr_ext_y and 10 * dy > dr_ext_y:
        ret += "Vertical Fresnel zone extension within [7,10]*pixel_height -> OK\n"
    else:
        ret += 'Vertical Fresnel zone extension NOT within [7,10]*pixel_height -> Check pixel width."\n'

    if Rx >= 3 * fwhm_x:
        ret += "Horizontal ROI > 3* FWHM(x) -> OK\n"
    else:
        ret += "Horizontal ROI !> 3* FWHM(x) -> Increase ROI width (x).\n"

    if Ry >= 3 * fwhm_y:
        ret += "Horizontal ROI > 3* FWHM(y) -> OK\n"
    else:
        ret += "Horizontal ROI !> 3* FWHM(y) -> Increase ROI height (y).\n"

    ret += "Focus sampling: FWHM > 10*px\n\n"

    ret += "END OF REPORT"

    return ret

def animate(wfr, qspace=False, logscale=False, delay = 10, outdir = None, fname = None):
    """ Generate an animated gif from the wavefront data. """
    intensity = wfr.get_intensity()

    # Get limits.
    xmin, xmax, ymax, ymin = wfr.get_limits()
    mx = intensity.max()
    mn = intensity.min()
    if logscale and mn <= 0.0:
        mn = intensity[np.where(intensity > 0.0)].min(),
        
    if fname is None:
        inp_filename = "wfr-animation"
    else:
        inp_filename = fname
        
    if outdir is None:
        os.mkdir("../../tmp/tmp")
        outdir = "../../tmp/"
        tmp_dir = "../../tmp/tmp"
    else:
        os.mkdir(outdir + "/tmp")
        tmp_dir = outdir + "/tmp"
        
    number_of_slices = intensity.shape[-1]
    # Setup a figure.

    for i in range(0,number_of_slices):

        #print("Processing slice #%d." % (i))

        # Plot profile as 2D colorcoded map.
        if logscale:
            plt.imshow(intensity[:,:,i], norm=mpl.colors.LogNorm(vmin=mn, vmax=mx), extent=[xmin*1.e6, xmax*1.e6, ymax*1.e6, ymin*1.e6], cmap="viridis")
        else:
            plt.imshow(intensity[:,:,i], norm=mpl.colors.Normalize(vmin=mn, vmax=mx), extent=[xmin*1.e6, xmax*1.e6, ymax*1.e6, ymin*1.e6], cmap="viridis")
        
        plt.xlabel("x ($\mu m$)")
        plt.ylabel("y ($\mu m$)")
        plt.savefig("%s/%s_%07d.png" % (tmp_dir, inp_filename, i) )
        plt.clf()

    os.system("convert -delay {} {}/*.png {}.gif".format(delay, outdir, inp_filename) )
    shutil.rmtree(tmp_dir)

def plotOnAxisPowerDensity(wfr, spectrum=False, outdir = None):
    """ Method to plot the on-axis power density.
    :param spectrum: Whether to plot the power density in energy domain (True) or time domain (False, default).
    :type spectrum: bool
    """
    """ Adapted from github:Samoylv/WPG/wpg/wpg_uti_wf.integral_intensity() """

    #print("\n Plotting on-axis power density.")
    # Setup new figure.
    plt.figure()

    # Switch to frequency (energy) domain if requested.
    if spectrum:
        srwl.SetRepresElecField(wfr._srwl_wf, 'f')
        intensity = wfr.get_intensity()

    # Get dimensions.
    mesh = wfr.params.Mesh
    dx = (mesh.xMax - mesh.xMin)/(mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin)/(mesh.ny - 1)

    # Get center pixel numbers.
    center_nx = int(mesh.nx/2)
    center_ny = int(mesh.ny/2)

    # Get time slices of intensity.
    intensity = wfr.get_intensity()

    # Get on-axis intensity.
    int0_00 = intensity[center_ny, center_nx, :]
    int0 = intensity.sum(axis=(0,1))*(dx*dy*1.e6) #  amplitude units sqrt(W/mm^2)
    int0max = int0.max()

    # Get meaningful slices.
    aw = [a[0] for a in numpy.argwhere(int0 > int0max*0.01)]
    if aw == []:
        raise RuntimeError("No significant intensities found.")

    dSlice = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices - 1)

    xs = numpy.arange(mesh.nSlices)*dSlice+ mesh.sliceMin
    xs_mf = numpy.arange(min(aw), max(aw))*dSlice + mesh.sliceMin

    # Plot.
    if(wfr.params.wDomain=='time'):
        plt.plot(xs*1e15,int0_00)
        plt.plot(xs_mf*1e15, int0_00[min(aw):max(aw)], 'ro')
        plt.title('On-Axis Power Density')
        plt.xlabel('time (fs)')
        plt.ylabel(r'Power density (W/mm${}^{2}$)')
    else: #frequency domain
        plt.plot(xs,int0_00)
        plt.plot(xs_mf, int0_00[min(aw):max(aw)], 'ro')
        plt.title('On-Axis Spectral Fluence')
        plt.xlabel('photon energy (eV)')
        plt.ylabel(r'fluence (J/eV/mm${}^{2}$)')

        # Switch back to time domain.
        srwl.SetRepresElecField(wfr._srwl_wf, 't')
        intensity = wfr.get_intensity()
    
    if outdir is not None:
        if spectrum == True:
            mode = 'spectrum'
        else:
            mode = 'time'
        plt.savefig(outdir + "/OnAxisPowerDensity_{}.png".format(mode))
        

def plotTotalPower(wfr, spectrum=False, outdir = None):
    """ Method to plot the total power.
    :param spectrum: Whether to plot the power density in energy domain (True) or time domain (False, default).
    :type spectrum: bool
    """

    """ Adapted from github:Samoylv/WPG/wpg/wpg_uti_wf.integral_intensity() """
    #print("\n Plotting total power.")
    # Setup new figure.
    plt.figure()

    # Switch to frequency (energy) domain if requested.
    if spectrum:
        #print("\n Switching to frequency domain.")
        srwl.SetRepresElecField(wfr._srwl_wf, 'f')
        intensity = wfr.get_intensity()
    else:
        intensity = wfr.get_intensity()
    # Get dimensions.
    mesh = wfr.params.Mesh
    dx = (mesh.xMax - mesh.xMin)/(mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin)/(mesh.ny - 1)

    # Get intensity by integrating over transverse dimensions.
    int0 = intensity.sum(axis=(0,1))

    # Scale to get unit W/mm^2
    int0 = int0*(dx*dy*1.e6) #  amplitude units sqrt(W/mm^2)
    int0max = int0.max()

    # Get center pixel numbers.
    center_nx = int(mesh.nx/2)
    center_ny = int(mesh.ny/2)

    # Get meaningful slices.
    aw = [a[0] for a in np.argwhere(int0 > int0max*0.01)]
    int0_mean = int0[min(aw):max(aw)]  # meaningful range of pulse
    dSlice = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices - 1)
    xs = np.arange(mesh.nSlices)*dSlice+ mesh.sliceMin
    xs_mf = np.arange(min(aw), max(aw))*dSlice + mesh.sliceMin
    if(wfr.params.wDomain=='time'):
        plt.plot(xs*1e15, int0) # time axis converted to fs.
        plt.plot(xs_mf*1e15, int0_mean, 'ro')
        plt.title('Power')
        plt.xlabel('time (fs)')
        plt.ylabel('Power (W)')
        dt = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices - 1)
        #print(('Pulse energy {:1.2g} J'.format(int0_mean.sum()*dt)))

    else: #frequency domain
        plt.plot(xs, int0)
        plt.plot(xs_mf, int0_mean, 'ro')
        plt.title('Spectral Energy')
        plt.xlabel('eV')
        plt.ylabel('J/eV')

        # Switch back to time domain.
        srwl.SetRepresElecField(wfr._srwl_wf, 't')
        wfr.intensity = wfr.get_intensity()
        
        if outdir is not None:
            if spectrum == True:
                mode = 'spectrum'
            else:
                mode = 'time'
            plt.savefig(outdir + "/TotalPower{}.png".format(mode))
            

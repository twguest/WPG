#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:03:50 2017

@author: gvanriessen
"""
# uncomment for Jupyter notebook
#%matplotlib notebook

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# Importing necessary modules:
import os
import sys

sys.path.insert(0, os.path.join("..", ".."))

import time
import copy
import numpy as np
import matplotlib.pyplot as plt
import itertools

import numexpr as ne


try:
    from wpg import srwlpy as srwl
except ImportError:
    import srwlpy as srwl  #  Hack for read the docs


# import SRW core functions
from wpg.srwlib import (
    SRWLOptD,
    SRWLOptA,
    SRWLOptC,
    SRWLOptT,
    SRWLOptL,
    SRWLOptMirEl,
    SRWLOptZP,
)

# import SRW helpers functions
# from wpg.useful_code.srwutils import AuxTransmAddSurfHeightProfileScaled


from wpg.useful_code.wfrutils import (
    calculate_fwhm_x,
    plot_wfront,
    calculate_fwhm_y,
    print_beamline,
    get_mesh,
    plot_1d,
    plot_2d,
)
from wpg.useful_code.wfrutils import propagate_wavefront

from wpg import Wavefront, Beamline
from wpg.generators import build_gauss_wavefront_xy  # Gaussian beam generator
from wpg.optical_elements import Empty, Use_PP
from wpg.wpg_uti_wf import *

from beautifultable import BeautifulTable


J2EV = 6.24150934e18


def timeit(method):
    # implements @timeit decorator:
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result

    return timed


def mkdir_p(path):
    """
    Create directory tree, if not exists (mkdir -p)

    :param path: Path to be created
    """
    if path == "":
        return
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def calculate_fwhm(wfr):
    """
    Calculate FWHM of the beam calculating number of point bigger then max/2 throuhgt center of the image

    :param wfr:  wavefront
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m]
    """
    intens = wfr.get_intensity(polarization="total").sum(axis=-1)

    mesh = wfr.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / mesh.nx
    dy = (mesh.yMax - mesh.yMin) / mesh.ny

    x_center = intens[intens.shape[0] // 2, :]
    fwhm_x = len(x_center[x_center > x_center.max() / 2]) * dx

    y_center = intens[:, intens.shape[1] // 2]
    fwhm_y = len(y_center[y_center > y_center.max() / 2]) * dy
    return {"fwhm_x": fwhm_x, "fwhm_y": fwhm_y}


def get_intensity_on_axis(wfr):
    """
    Calculate intensity (spectrum in frequency domain) along z-axis (x=y=0)

    :param wfr:  wavefront
    :return: [z,s0] in [a.u.] if frequency domain
    """

    wf_intensity = wfr.get_intensity(polarization="horizontal")
    mesh = wfr.params.Mesh
    zmin = mesh.sliceMin
    zmax = mesh.sliceMax
    sz = numpy.zeros((mesh.nSlices, 2), dtype="float64")
    sz[:, 0] = numpy.linspace(zmin, zmax, mesh.nSlices)
    sz[:, 1] = wf_intensity[mesh.nx / 2, mesh.ny / 2, :] / wf_intensity.max()

    return sz


def objOpt(obj):
    # wrapper function for  cleaning up code  a little bit - replace with custom class later.
    from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file

    opt = srwl_opt_setup_transm_from_file(
        obj.file_path,
        obj.resolution,
        obj.thickness,
        obj.delta,
        obj.atten_len,
        xc=obj.xc if obj.contains("xc") else 0.0,
        yc=obj.yc if obj.contains("yc") else 0.0,
        area=obj.area if obj.contains("area") else None,
        rotate_angle=obj.rotate_angle if obj.contains("rotate_angle") else None,
        rotate_reshape=obj.rotate_reshape if obj.contains("rotate_reshape") else None,
        cutoff_background_noise=obj.cutoff_background_noise
        if obj.contains("cutoff_background_Noise")
        else 0,
        background_color=obj.background_color
        if obj.contains("background_color")
        else 0,
        tile=obj.tile if obj.contains("tile") else None,
        shift_x=obj.shift_x if obj.contains("shift_x") else None,
        shift_y=obj.shift_y if obj.contains("shift_y") else None,
        invert=obj.invert if obj.contains("invert") else None,
        is_save_images=obj.is_save_images if obj.contains("is_save_images") else None,
        prefix=obj.prefix if obj.contains("prefix") else None,
        output_image_format=obj.output_image_format
        if obj.contains("output_image_format")
        else None,
    )

    return opt


def snakeScan(centre=(0, 0), range=(10e-6, 10e-6), spacing=(1.0e-6, 1.0e-6)):

    from pylab import arange, concatenate, repeat

    if len(spacing) == 1:
        spacing = (spacing, spacing)

    # define some grids
    xgrid = arange(20, 31)
    ygrid = arange(10, 16)

    xscan = concatenate([xgrid[:: (-1) ** i] for i in range(len(ygrid))])
    yscan = repeat(ygrid, len(xgrid))
    return list(zip(a, b))


# copied from Pynx
def spiral_archimedes(a, n):
    """" Creates np points spiral of step a, with a between successive points
    on the spiral. Returns the x,y coordinates of the spiral points.

    This is an Archimedes spiral. the equation is:
      r=(a/2*pi)*theta
      the stepsize (radial distance between successive passes) is a
      the curved absciss is: s(theta)=(a/2*pi)*integral[t=0->theta](sqrt(1*t**2))dt
    """
    vr, vt = [0], [0]
    t = np.pi
    while len(vr) < n:
        vt.append(t)
        vr.append(a * t / (2 * np.pi))
        t += 2 * np.pi / np.sqrt(1 + t ** 2)
    vt, vr = np.array(vt), np.array(vr)
    return vr * np.cos(vt), vr * np.sin(vt)


# copied from Pynx
def spiral_fermat(dmax, n):
    """"
    Creates a Fermat spiral with n points distributed in a circular area with
    diameter<= dmax. Returns the x,y coordinates of the spiral points. The average
    distance between points can be roughly estimated as 0.5*dmax/(sqrt(n/pi))

    http://en.wikipedia.org/wiki/Fermat%27s_spiral
    """
    c = 0.5 * dmax / np.sqrt(n)
    vr, vt = [], []
    t = 0.4
    goldenAngle = np.pi * (3 - np.sqrt(5))
    while t < n:
        vr.append(c * np.sqrt(t))
        vt.append(t * goldenAngle)
        t += 1
    vt, vr = np.array(vt), np.array(vr)
    return vr * np.cos(vt), vr * np.sin(vt)


class Struct:
    "A structure that can have any fields defined."

    def __init__(self, **entries):
        self.__dict__.update(entries)

    def contains(self, attrName):
        return hasattr(self, attrName)


def get_transmissionFunctionCplx(tr):

    """
    Extract complex transmission function of  transmission object

    :param transmission: SRWLOptT struct, see srwlib.h
    :return: mag, phase tuple map of transmission object
    """

    mesh = tr.mesh
    nx = mesh.nx
    ny = mesh.ny
    phase = np.array(tr.arTr[1::2]).reshape((ny, nx))
    amplitude = np.array(tr.arTr[::2]).reshape((ny, nx))

    cplx = ne.evaluate("complex(phase, amplitude)")

    return cplx


def add_transmissionFunctions(trA, trB):

    return ne.evaluate("trA + trB")


def show_transmissionFunctionCplx(trCplx):
    from husl import complex_to_rgb

    rgb = complex_to_rgb(
        z=trCplx,
        amin=None,
        amax=None,
        mode="special",
        phstart=0.0,
        sat=1.0,
        as_image=True,
    )

    # make colour wheel:
    colorWheel = complex_to_rgb(
        z=None,
        amin=None,
        amax=None,
        mode="special",
        phstart=0.0,
        sat=1.0,
        as_image=True,
    )

    plt.figure(figsize=(10, 7))
    plt.subplot(121)
    plt.imshow(rgb)
    plt.title("Transmission Function")
    plt.xlabel("mm")
    plt.ylabel("mm")
    plt.subplot(122)
    plt.imshow(colorWheel)
    plt.title("Colorwheel (TODO: proper labelling")
    plt.xlabel("")
    plt.ylabel("")
    plt.show()


# function to allow more convenient use of default propgation parameters)
# all instances of this should be replaced with class Use_PP(object) in opticalElements.py:
def propagationParameters(
    AutoResizeBeforeProp=0,  # 0
    AutoResizeAfterProp=0,  # 1
    RelPrecAutoresize=1.0,  # 2
    SemiAnalyt=0,  # 3
    ResizeFourierSide=0,  # 4
    RangeX=1.0,  # 5
    ResolutionX=1.0,  # 6
    RangeY=1.0,  # 7
    ResolutionY=1.0,  # 8
    ShiftType=0,  # 9
    ShiftX=0,  # 10
    ShiftY=0,
):  # 11

    # propagation parameter list & description:
    # [0]:  Auto-Resize (1) or not (0) Before propagation
    # [1]:  Auto-Resize (1) or not (0) After propagation
    # [2]:  Relative Precision for propagation with Auto-Resizing (1. is nominal)
    # [3]:  Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation
    # [4]:  Do any Resizing on Fourier side, using FFT, (1) or not (0)
    # [5]:  Horizontal Range modification factor at Resizing (1. means no modification)
    # [6]:  Horizontal Resolution modification factor at Resizing
    # [7]:  Vertical Range modification factor at Resizing
    # [8]:  Vertical Resolution modification factor at Resizing
    # [9]:  Type of wavefront Shift before Resizing (not yet implemented)
    # [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    # [11]: New Vertical wavefront Center position after Shift (not yet implemented)

    pp = [
        AutoResizeBeforeProp,
        AutoResizeAfterProp,
        RelPrecAutoresize,
        SemiAnalyt,
        ResizeFourierSide,
        RangeX,
        ResolutionX,
        RangeY,
        ResolutionY,
        ShiftType,
        ShiftX,
        ShiftY,
    ]

    return pp


def pixelScale(lengthInMeters, metersPerPixel):
    return lengthInMeters / metersPerPixel


def joinContainers(opticalContainerList):
    # concatenates list of SRWLOptC objects
    # note that the beamline.append method could can be used instead of this
    oe = []
    pp = []
    for ocl in opticalContainerList:
        if ocl != None:
            oe.append(ocl.getOpticalElements())
            pp.append(ocl.getPropagationParameters())

    container = SRWLOptC(_arOpt=oe, _arProp=pp)

    return container


def visVol(volData):
    """ This example demonstrates rendering a color volume.
    We demonstrate two renderers capable of rendering color data:
    the colormip and coloriso renderer.
    """

    import visvis as vv

    app = vv.use()

    # Load volume
    vol = volData

    # set labels
    vv.xlabel("x axis")
    vv.ylabel("y axis")
    vv.zlabel("z axis")

    #    # Create figure and make subplots with different renderers
    #    vv.figure(1); vv.clf()
    #    RS = ['mip', 'iso', 'edgeray', 'ray', 'litray']
    #    a0 = None
    #    tt = []
    #    for i in range(5):
    #        a = vv.subplot(3,2,i+2)
    #        t = vv.volshow(vol)
    #        vv.title('Renderstyle ' + RS[i])
    #        t.colormap = vv.CM_HOT
    #        t.renderStyle = RS[i]
    #        t.isoThreshold = 200  # Only used in iso render style
    #        tt.append(t)
    #        if a0 is None:
    #            a0 = a
    #        else:
    #            a.camera = a0.camera
    t = vv.volshow(vol, renderStyle="edgeray")
    t.colormap = vv.CM_HOT

    # Get axes and set camera to orthographic mode (with a field of view of 70)
    a = vv.gca()
    a.camera.fov = 45

    # Create colormap editor wibject.
    vv.ColormapEditor(a)

    # Create colormap editor in first axes
    # cme = vv.ColormapEditor(vv.gcf(), t)

    # Run app
    # app.Create()
    app.Run()


def calc_sampling(zoom, mf):
    """
    This function calculates sampling.
    :param zoom: range zoom
    :param mf: modification factor for step, i.e. dx1=mf*dx0
    :return: sampling.
    """
    sampling = zoom / mf
    print("zoom:{:.1f}; mod_factor:{:.1f}; sampling:{:.1f}".format(zoom, mf, sampling))
    return sampling


def _resample(wf, axis, data, x0, x1):
    if axis.lower() == "x":
        y = data[data.shape[0] / 2, :]
        x = np.linspace(wf.params.Mesh.xMin, wf.params.Mesh.xMax, y.shape[0])
    elif axis.lower() == "y":
        y = data[:, data.shape[1] / 2]
        x = np.linspace(wf.params.Mesh.yMin, wf.params.Mesh.yMax, y.shape[0])
    else:
        raise ValueError('Wrong axis {}, should be "x" or "y"'.format(axis))

    if not x0 is None:
        xmin = x0
    else:
        xmin = x[0]

    if not x1 is None:
        xmax = x1
    else:
        xmax = x[-1]

    x1 = np.linspace(xmin, xmax, len(y))
    y1 = np.interp(x1, x, y)
    return x1, y1


def qParameter(PhotonEnergy, Waist, RadiusCurvature):
    """Computing complex q parameter
    see test_SRWLIB_Example15.py"""
    Lam = 1.24e-6 * PhotonEnergy
    qp = (1.0 + 0j) / complex(1 / RadiusCurvature, -Lam / 3.1415 / Waist ** 2)
    return qp, Lam


def intensity_cut(wf, axis, polarization, x0=None, x1=None):

    if polarization.lower() == "v" or polarization.lower() == "vertical":
        pol = "vertical"
    elif polarization.lower() == "h" or polarization.lower() == "horizontal":
        pol = "horizontal"
    elif polarization.lower() == "t" or polarization.lower() == "total":
        pol = "total"
    else:
        raise ValueError(
            'Wrong polarization {}, should be "v" or "vertical"'
            + ' or "h" or "horizontal" or "t" or "total"'.format(polarization)
        )

    data = wf.get_intensity(slice_number=0, polarization=pol)
    return _resample(wf, axis, data, x0, x1)


def phase_cut(wf, axis, polarization, x0=None, x1=None):

    if polarization.lower() == "v" or polarization.lower() == "vertical":
        pol = "vertical"
    elif polarization.lower() == "h" or polarization.lower() == "horizontal":
        pol = "horizontal"
    else:
        raise ValueError(
            'Wrong polarization {}, should be "v" or "vertical" or "h" or "horizontal"'.format(
                polarization
            )
        )

    data = wf.get_phase(slice_number=0, polarization=pol)
    return _resample(wf, axis, data, x0, x1)


def defineEFM(orient, p, q, thetaEFM, theta0, lengthEFM):
    """
    A wrapper to a SRWL function SRWLOptMirEl() for defining a plane elliptical focusing mirror propagator

    :param Orient:    mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param p:  the distance to two ellipsis centers
    :param q:  the distance to two ellipsis centers
    :param thetaEFM:  the design incidence angle in the center of the mirror
    :param theta0:    the "real" incidence angle in the center of the mirror
    :param lengthEFM: mirror length, [m]
    :return: the struct opEFM
    """
    if orient == "x":  # horizontal plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(
            _p=p,
            _q=q,
            _ang_graz=thetaEFM,
            _r_sag=1.0e40,
            _size_tang=lengthEFM,
            _nvx=np.cos(theta0),
            _nvy=0,
            _nvz=-np.sin(theta0),
            _tvx=-np.sin(theta0),
            _tvy=0,
            _x=0,
            _y=0,
            _treat_in_out=1,
        )
    elif orient == "y":  # vertical plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(
            _p=p,
            _q=q,
            _ang_graz=thetaEFM,
            _r_sag=1.0e40,
            _size_tang=lengthEFM,
            _nvx=0,
            _nvy=np.cos(theta0),
            _nvz=-np.sin(theta0),
            _tvx=0,
            _tvy=-np.sin(theta0),
            _x=0,
            _y=0,
            _treat_in_out=1,
        )
    else:
        raise TypeError('orient should be "x" or "y"')
    return opEFM


def meshDim(wf):
    """
    return wavefront mesh dimensions
    """
    wf_mesh = wf.params.Mesh
    return wf_mesh.nx, wf_mesh.ny


def calculate_source_fwhm(ekev, theta_fwhm):
    """
    Calculate source size from photon energy and FWHM angular divergence

    :param evev: Energy in keV
    :param theta_fwhm: theta_fwhm [units?]
    """
    wl = 12.39e-10 / ekev
    k = 2 * np.sqrt(2 * np.log(2))
    theta_sigma = theta_fwhm / k
    sigma0 = wl / (2 * np.pi * theta_sigma)
    return sigma0 * k


def calculate_theta_fwhm_cdr(ekev, qnC):
    """
    Calculate angular divergence using formula from XFEL CDR2011

    :param ekev: Energy in keV
    :param qnC: e-bunch charge, [nC]
    :return: theta_fwhm [units?]
    """
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC)) * 1e-6 / ekev ** 0.85
    return theta_fwhm


def modes2D(m, N):
    """
    return N pairs of Transverse Gauss-Hermite Mode Order pairs with order up to m.

    """
    A = list()
    for i in range(1, m + 1):
        A.append(list(itertools.product([0, i], repeat=2)))

    # combine
    A = list(itertools.chain.from_iterable(A))

    # remove duplicates
    temp = []
    for a, b in A:
        if (a, b) not in temp:  # to check for the duplicate tuples
            temp.append((a, b))

    return temp[0:N]


def eigenvaluePartCoh(p_I, p_mu, n):
    """
    GVR

    return eigenvalue normalised to eigenvalue of fundamental

    definitions follow Starikov and Wolf, 1982:


    p_mu and p_I are the rms widths of the degree of coherence and of the intensity of the source.

    beta =  p_mu/p_I is a measure of the "degree of global coherence" of the source.

    When beta >> 1, the source is effectively spatially coherent in the global sense and is then
    found to be well represented by a single mode.
    When beta << 1, the source is effectively spatially incoherent in the global
    sense, and the number of modes needed to describe its behavior is of the order of beta/3

    """

    a = 1 / (2 * p_I ** 2)
    b = 1 / (2 * p_mu ** 2)
    c = (a ** 2 + 2 * a * b) ** (1 / 2)
    l_0 = 1.0

    l_n = l_0 * (b / (a + b + c)) ** n

    return l_n


def plotEigenValues(e, _threshold):

    plt.plot(e, "bo")
    plt.plot(e[e > _threshold], "r+")
    plt.title("Eigenvalues")
    plt.xlabel("Mode")
    plt.ylabel("$\lambda/\lambda_0$")
    plt.legend(loc=2)
    plt.show()


@timeit
def plotWavefront(
    wf, title, slice_numbers=False, cuts=False, interactive=False, phase=False
):
    # draw wavefront with common functions
    # if phase = True, plot phase only

    print("Displaying wavefront...")

    if slice_numbers is None:
        slice_numbers = range(wf_intensity.shape[-1])

    if isinstance(slice_numbers, int):
        slice_numbers = [
            slice_numbers,
        ]

    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wf)
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)
    print("stepX, stepY [um]:", dx * 1e6, dy * 1e6, "\n")

    if phase == True:
        A = wf.get_phase(slice_number=0, polarization="horizontal")
        label = "Phase (rad.)"
    else:
        ii = wf.get_intensity(slice_number=0, polarization="horizontal")
        ii = ii * wf.params.photonEnergy / J2EV  # *1e3
        imax = numpy.max(ii)

        if wf.params.wEFieldUnit != "arbitrary":
            print(
                "Total power (integrated over full range): %g [GW]"
                % (ii.sum(axis=0).sum(axis=0) * dx * dy * 1e6 * 1e-9)
            )
            print(
                "Peak power calculated using FWHM:         %g [GW]"
                % (
                    imax
                    * 1e-9
                    * 1e6
                    * 2
                    * numpy.pi
                    * (calculate_fwhm_x(wf) / 2.35)
                    * (calculate_fwhm_y(wf) / 2.35)
                )
            )
            print("Max irradiance: %g [GW/mm^2]" % (imax * 1e-9))
            label = "Irradiance (W/$mm^2$)"
        else:
            ii = ii / imax
            label = "Irradiance (a.u.)"

        A = ii

    [x1, x2, y1, y2] = wf.get_limits()

    if interactive == True:
        import mpld3
        from mpld3 import plugins
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        im = ax.imshow(ii, extent=[x1 * 1e3, x2 * 1e3, y1 * 1e3, y2 * 1e3])
        fig.colorbar(im, ax=ax)

        fig = plt.figure()
        ax.set_title(title, size=20)
        plugins.connect(fig, plugins.MousePosition(fontsize=14))
        mpld3.show()
        # mpld3.display(fig)

    else:
        pylab.figure(figsize=(10, 6))
        pylab.imshow(A, extent=[x1 * 1e3, x2 * 1e3, y1 * 1e3, y2 * 1e3])
        pylab.set_cmap("hot")
        pylab.axis("tight")
        pylab.colorbar(orientation="horizontal")
        pylab.xlabel("x (mm)")
        pylab.ylabel("y (mm)")
        pylab.axes().set_aspect(0.5)

        pylab.title(title)
        pylab.show()


#
#
#    if wf.params.Mesh.nSlices==1:
#        dt = 1  # added because alternative results in in div by zero where only 1 slice exists
#    else:
#        dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / (wf.params.Mesh.nSlices - 1)
#
#    for sn in slice_numbers:
#        data = wf_intensity[:, :, sn]
#        data = data*dt
#        phase = wf_phase[:, :, sn]
#
#        plt.subplot(1,2,1)
#        plt.imshow(data)
#        plt.subplot(1,2,2)
#        plt.imshow(phase)
#        plt.show()
#
#    #draw wavefront with cuts
#    if cuts:
#        plot_wfront2(wf, title_fig=title,
#                    isHlog=False, isVlog=False,
#                    orient='x', onePlot=True,
#                    i_x_min=None, i_y_min=None)
#
#        plt.set_cmap('jet') #set color map, 'bone', 'hot', 'jet', etc
#        plt.show()


def plot_wfront2(
    mwf,
    title_fig,
    isHlog,
    isVlog,
    orient,
    onePlot,
    i_x_min=None,
    i_y_min=None,
    bPlotPha=None,
    saveDir=None,
):
    """ Adapted from plot_wfront in wfrutils


        Plot 2D wavefront (a slice).

        :param mwf: 2D wavefront structure
        :param title_fig: Figure title
        :param isHlog: if True, plot the horizontal cut in logarithmic scale
        :param isVlog: if True, plot the vertical cut in logarithmic scale
        :param i_x_min: Intensity threshold for horizontral cut, i.e. x-axis  limits are [min(where(i_x<i_x_min):max(where(i_x<i_x_min)]
        :param i_y_min: Intensity threshold for vertical cut,
        :param orient: 'x' for returning horizontal cut, 'y' for vertical cut
        :param onePlot: if True, put intensity map and plot of cuts on one plot, as  subplots
        :param bPlotPha: if True, plot the cuts of WF phase
        :return: 2-column array containing horizontal or vertical cut data in dependence of 'orient' parameter
        """
    if isHlog:
        print("FWHMx[um]:", calculate_fwhm_x(mwf) * 1e6)
    else:
        print("FWHMx [mm]:", calculate_fwhm_x(mwf) * 1e3)
    if isVlog:
        print("FWHMy [um]:", calculate_fwhm_y(mwf) * 1e6)
    else:
        print("FWHMy [mm]:", calculate_fwhm_y(mwf) * 1e3)
    [xc, yc] = calculate_peak_pos(mwf)
    print("Coordinates of center, [mm]:", xc * 1e3, yc * 1e3)
    ii = mwf.get_intensity(slice_number=0, polarization="horizontal")
    # [LS14-06-02]
    # for 2D Gaussian the intrincic SRW GsnBeam wave field units Nph/mm^2/0.1%BW
    # to get fluence W/mm^2
    # (Note: coherence time for Gaussian beam duration should be specified):
    ii = ii * mwf.params.photonEnergy / J2EV  # *1e3
    imax = numpy.max(ii)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    ph = mwf.get_phase(slice_number=0, polarization="horizontal")
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)
    print("stepX, stepY [um]:", dx * 1e6, dy * 1e6, "\n")
    xa = numpy.linspace(xmin, xmax, nx)
    ya = numpy.linspace(ymin, ymax, ny)

    if mwf.params.wEFieldUnit != "arbitrary":
        print(
            "Total power (integrated over full range): %g [GW]"
            % (ii.sum(axis=0).sum(axis=0) * dx * dy * 1e6 * 1e-9)
        )
        print(
            "Peak power calculated using FWHM:         %g [GW]"
            % (
                imax
                * 1e-9
                * 1e6
                * 2
                * numpy.pi
                * (calculate_fwhm_x(mwf) / 2.35)
                * (calculate_fwhm_y(mwf) / 2.35)
            )
        )
        print("Max irradiance: %g [GW/mm^2]" % (imax * 1e-9))
        label4irradiance = "Irradiance (W/$mm^2$)"
    else:
        ii = ii / imax
        label4irradiance = "Irradiance (a.u.)"

    pylab.figure(figsize=(21, 6))
    if onePlot:
        pylab.subplot(131)
    [x1, x2, y1, y2] = mwf.get_limits()
    pylab.imshow(ii, extent=[x1 * 1e3, x2 * 1e3, y1 * 1e3, y2 * 1e3])
    pylab.set_cmap("bone")
    # pylab.set_cmap('hot')
    pylab.axis("tight")
    # pylab.colorbar(orientation='horizontal')
    pylab.xlabel("x (mm)")
    pylab.ylabel("y (mm)")
    pylab.title(title_fig)

    irr_y = ii[:, numpy.max(numpy.where(xa == xc))]
    irr_x = ii[numpy.max(numpy.where(ya == yc)), :]
    pha_y = ph[:, numpy.max(numpy.where(xa == xc))]
    pha_x = ph[numpy.max(numpy.where(ya == yc)), :]

    if i_y_min is None:
        i_y_min = numpy.min(irr_x)

    if i_x_min is None:
        i_x_min = numpy.min(irr_x)

    if onePlot:
        pylab.subplot(132)
    else:
        pylab.figure()
    if isVlog and numpy.max(irr_y) > 0:
        # ya = ya*1e6
        pylab.semilogy(ya * 1e6, irr_y, "-vk")
        pylab.xlabel("(um)")
        pylab.xlim(
            numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e6,
            numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e6,
        )
    else:
        # ya = ya*1e3
        pylab.plot(ya * 1e3, irr_y)
        pylab.xlabel("y (mm)")
        pylab.xlim(
            numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3,
            numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3,
        )
    pylab.ylim(0, numpy.max(ii) * 1.1)
    pylab.ylabel(label4irradiance)
    pylab.title("Vertical cut,  xc = " + str(int(xc * 1e6)) + " um")
    pylab.grid(True)
    if onePlot:
        pylab.subplot(133)
    else:
        pylab.figure()
    if isHlog and numpy.max(irr_x) > 0:
        # xa = xa*1e6
        pylab.semilogy(xa * 1e6, irr_x, "-vr")
        pylab.xlabel("x, (um)")
        pylab.xlim(
            numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e6,
            numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e6,
        )
    else:
        # xa = xa*1e3
        pylab.plot(xa * 1e3, irr_x)
        pylab.xlabel("x (mm)")
        pylab.xlim(
            numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3,
            numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3,
        )
    pylab.ylim(0, numpy.max(ii) * 1.1)
    pylab.ylabel(label4irradiance)
    pylab.title("Horizontal cut, yc = " + str(int(yc * 1e6)) + " um")
    pylab.grid(True)

    if saveDir is not None:
        epsname = "%s/%s.eps" % (saveDir, title_fig.split("at ")[1].split(" m")[0])
        pylab.savefig(epsname)
        # pylab.close(epsfig)

    if bPlotPha:
        pylab.figure()
        pylab.plot(ya * 1e3, pha_y, "-ok")
        pylab.xlim(
            numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3,
            numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3,
        )
        pylab.ylim(-numpy.pi, numpy.pi)
        pylab.xlabel("y (mm)")
        pylab.title("phase, vertical cut, x=0")
        pylab.grid(True)

        pylab.figure()
        pylab.plot(xa * 1e3, pha_x, "-or")
        pylab.xlim(
            numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3,
            numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3,
        )
        pylab.ylim(-numpy.pi, numpy.pi)
        pylab.xlabel("x (mm)")
        pylab.title("phase, horizontal cut, y=0")
        pylab.grid(True)

    if orient == "x":
        dd = numpy.zeros(shape=(nx, 2), dtype=float)
        dd[:, 0] = xa
        # for idx in range(nx): dd[idx, 1] = sum(ii[:, idx])
        dd[:, 1] = irr_x
    if orient == "y":
        dd = numpy.zeros(shape=(ny, 2), dtype=float)
        dd[:, 0] = ya
        # for idx in range(ny): dd[idx, 1] = sum(ii[idx, :])
        dd[:, 1] = irr_y
    return dd


def ZPf(D, drn, E, _n=1):

    """
    D is outer diamter [m]
    drn is width  of outer zone [m]
    E is energy [eV]
    _n is diffraction order (default 1])
    """

    w = 1.2398 / E  # use swrlib swrl_uti_ph_en_conv in future!
    f = (D * (drn * 1.0e6)) / (_n * w)
    return f


def writeCSV(wf, fn, _comment=""):

    np.savetxt(fn, wf, delimiter=",", comments=_comment)


def animateWF(wf, _filename=""):
    """
    wf should be a wf with multiple slices... For now (Testing) it is a list of image arrays
    if filename given a movie file will be saved
    """

    import matplotlib.animation as animation

    fig = plt.figure()

    ims = []

    for i in range(0, len(wf) - 1):
        print(i)

        im = plt.imshow(wf[i], animated=True)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=500)
    plt.show()

    if _filename != "":
        ani.save(_filename)


def constructWavefield(
    ekeV,
    qnC,
    z1,
    range_xy,
    strOutputDataFolder=None,
    dimension=2,
    beta=0.7,  # 'tunable' parameter — sets global degree of coherence
    threshold=0.8,  # ignore modes with eigenvalue < this value
    npoints=512,
    nslices=10,
    display=True,
):
    """
    Define Transverse Optical Modes
    """
    k = 2 * np.sqrt(2 * np.log(2))
    wlambda = 12.4 * 1e-10 / ekeV  # wavelength
    theta_fwhm = calculate_theta_fwhm_cdr(ekeV, qnC)
    sigX = 12.4e-10 * k / (ekeV * 4 * np.pi * theta_fwhm)

    # coherence width:
    widthCoherence = beta * sigX

    # get first n eigenvalues for transverse modes
    n = 99
    e0 = eigenvaluePartCoh(sigX, widthCoherence, range(0, n))

    # keep only modes for which eigenvalue > threshold
    e = e0[e0 > threshold]

    if display:
        plotEigenValues(e0, threshold)

    # generate mode indices
    modes = modes2D(9, N=len(e))

    dimension = 2  # value should be 3 for 3D wavefront, 2 for 2D wavefront
    wf = []
    for mx, my in modes:

        # define unique filename for storing results
        ip = np.floor(ekeV)
        frac = np.floor((ekeV - ip) * 1e3)

        # build initial gaussian wavefront
        if dimension == 2:
            wfr0 = build_gauss_wavefront_xy(
                nx=npoints,
                ny=npoints,
                ekev=ekeV,
                xMin=-range_xy / 2,
                xMax=range_xy / 2,
                yMin=-range_xy / 2,
                yMax=range_xy / 2,
                sigX=sigX,
                sigY=sigX,
                d2waist=z1,
                _mx=mx,
                _my=my,
            )

        else:
            # build initial 3d gaussian beam
            tau = 1
            # not sure if this parameter is even used - check meaning.
            wfr0 = build_gauss_wavefront(
                nx=npoints,
                ny=npoints,
                nz=nslices,
                ekev=ekev,
                xMin=-range_xy / 2,
                xMax=range_xy / 2,
                yMin=-range_xy / 2,
                yMax=range_xy / 2,
                tau=tau,
                sigX=sigX,
                sigY=sigX,
                d2waist=z1,
                _mx=mx,
                _my=my,
            )
            if display == True:
                print(
                    "dy {:.1f} um".format(
                        (mwf.params.Mesh.yMax - mwf.params.Mesh.yMin)
                        * 1e6
                        / (mwf.params.Mesh.ny - 1.0)
                    )
                )
                print(
                    "dx {:.1f} um".format(
                        (mwf.params.Mesh.xMax - mwf.params.Mesh.xMin)
                        * 1e6
                        / (mwf.params.Mesh.nx - 1.0)
                    )
                )
                plot_t_wf(mwf)
                look_at_q_space(mwf)

        # init WPG Wavefront helper class
        mwf = Wavefront(wfr0)

        # store wavefront to HDF5 file
        if strOutputDataFolder:
            fname0 = (
                "g"
                + str(int(ip))
                + "_"
                + str(int(frac))
                + "kev"
                + "_tm"
                + str(mx)
                + str(my)
            )
            ifname = os.path.join(strOutputDataFolder, fname0 + ".h5")
            mwf.store_hdf5(ifname)
            print("Saved wavefront to HDF5 file: {}".format(ifname))
        else:
            ifname = None

        wf.append([mx, my, ifname, mwf])

        # plotWavefront(mwf, 'at '+str(z1)+' m')
        # look_at_q_space(mwf)

        fwhm_x = calculate_fwhm_x(mwf)
        print(
            "FWHMx [mm], theta_fwhm [urad]: {}, {}".format(
                fwhm_x * 1e3, fwhm_x / z1 * 1e6
            )
        )

        # show_slices_hsv(mwf, slice_numbers=None, pretitle='SLICETYSLICE')

        return wf, modes, e


def writeIntensity(
    wavefield, label, path="./out", polarization="horizontal", imgType="tif"
):

    from slugify import slugify
    import imageio
    import os

    # imageio.imwrite(os.path.join(path, slugify(label)+'.'+imgType),
    #                wavefield.get_intensity(polarization=polarization))

    imageio.imwrite(
        os.path.join(path, label + "." + imgType),
        wavefield.get_intensity(polarization=polarization),
    )

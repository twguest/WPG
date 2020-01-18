# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

__author__ = 'A. Buzmakov'

import warnings
from wpg import srwlib

try:
    from wpg import srwlpy
except ImportError:
    import srwlpy  #  Hack for read the docs

def build_gaussian(GsnBm, nx, ny, xMin, xMax, yMin, yMax):
    """
    Build 2D Gaussian beam.
    
    :param nx: Number of point along x-axis
    :param ny: Number of point along y-axis
    :param nz: Number of point along z-axis (slices)
    :param ekev: Energy in kEv
    :param xMin: Initial Horizontal Position [m]
    :param xMax: Final Horizontal Position [m]
    :param yMin: Initial Vertical Position [m]
    :param yMax: Final Vertical Position [m]
    :param sigX: Horiz. RMS size at Waist [m]
    :param sigY:  Vert. RMS size at Waist [m]
    :param d2waist: distance to Gaussian waist
    :param xoff:    Horizonal Coordinate of Gaussian Beam Center at Waist [m]
    :param yoff:    Vertical  Coordinate of Gaussian Beam Center at Waist [m]
    :param tiltX:   Average Angle of Gaussian Beam at Waist in Horizontal plane [rad] 
    :param tiltY:   Average Angle of Gaussian Beam at Waist in Vertical plane [rad]
    :param pulseEn:  Energy per Pulse [J]
    :param pulseTau: Coherence time [s] to get proper BW
    :param _mx: transverse Gauss-Hermite mode order in horizontal direction
    :param _my: transverse Gauss-Hermite mode order in vertical direction
    :return: wpg.Wavefront structure

    """
    wfr = srwlib.SRWLWfr()  # Initial Electric Field Wavefront
    wfr.allocate(1, nx, ny)
     # Numbers of points vs Photon Energy (1), Horizontal and
     # Vertical Positions (dummy)
    wfr.mesh.eStart = GsnBm.avgPhotEn  # Initial Photon Energy [eV]
    wfr.mesh.eFin = GsnBm.avgPhotEn  # Final Photon Energy [eV]
    wfr.avgPhotEn = (wfr.mesh.eStart + wfr.mesh.eFin) / 2
    wfr.mesh.zStart = 0
    wfr.mesh.xStart = xMin  # Initial Horizontal Position [m]
    wfr.mesh.xFin = xMax  # Final Horizontal Position [m]
    wfr.mesh.yStart = yMin  # Initial Vertical Position [m]
    wfr.mesh.yFin = yMax  # Final Vertical Position [m]

    #wfr.presFT = 1  # Defining Initial Wavefront in Time Domain
    wfr.presFT = 0  # Defining Initial Wavefront in Freq Domain

    # Some information about the source in the Wavefront structure
    wfr.partBeam.partStatMom1.x = GsnBm.x
    wfr.partBeam.partStatMom1.y = GsnBm.y
    wfr.partBeam.partStatMom1.z = GsnBm.z
    wfr.partBeam.partStatMom1.xp = GsnBm.xp
    wfr.partBeam.partStatMom1.yp = GsnBm.yp


    sampFactNxNyForProp = -1  # sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [sampFactNxNyForProp]
    srwlpy.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
    return wfr

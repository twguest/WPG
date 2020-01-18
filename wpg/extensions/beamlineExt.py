#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:21:08 2019

@author: gvanriessen
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""
This module contains wrapper for SRWLOptC (optical container) and WPG beamline

"""

import wpg.srwlib as srwlib
from wpg.srwlib import srwl
from wpg.utils import srw_obj2str
from wpg.optical_elements import Screen, Empty, WPGOpticalElement
from wpg.wpg_uti_wf import plot_intensity_map, averaged_intensity,calculate_fwhm, get_intensity_on_axis
from wpg.wpg_uti_oe import show_transmission 
import os.path

from wpg.beamline import Beamline
from extensions.gvrutils import plotWavefront, show_transmissionFunctionCplx, J2EV

from wpg.srwlib import SRWLOptD as Drift

from extensions.husl import complex_to_rgb
from wpg.useful_code.wfrutils import get_mesh
import imageio
import numpy as np

class bl(Beamline):
    """
    Set of optical elements and propagation parameters.
    """

    def __init__(self, srwl_beamline=None, description=None):
         super(bl, self).__init__(srwl_beamline=srwl_beamline)
                
         self.description=description       


    def __str__(self):
        """
        String representaion of beamline (used with print function).

        :return: string
        """        
        return self.description + ': ' + super(bl, self).__str__()
  
  
    def propagate(self, 
                  wfr,  
                  plot=False, 
                  describe=True,
                  outFilePath=None, 
                  propagationOverride=None):
        """
        Propagate wavefront through beamline.

        :param wfr: Input wavefront (will be rewritten after propagation)
        :type wfr: wpg.wavefront.Wavefront
        """
        
        print ('Propagating through ' + self.description)
        super(bl, self).propagate(wfr)
        
        if (plot == True) or (describe == True) or (outFilePath is not None):
          
            # Create a temporary screen 
            screen = oeScreen(filename=outFilePath,
                              writeIntensity=(outFilePath!=None),
                              showIntensity=plot,
                              description=self.description + ' [screen]',                              
                              writePhase = False,
                              showPhase = False,
                              calculateMetrics = False)
            metrics = screen.propagate(wfr)
        else:
            metrics = None
            
        
        #print ('Intensity on axis = %e ' % (get_intensity_on_axis(wfr)))
        #print ('FWHM = %e ' % (calculate_fwhm(wfr)))
        return metrics
 
#class oeDrift(Drift):
#    
#    """Optical Element: Drift Space"""
#    
#    def __init__(self,length=0, description=None):
#        """
#        :param _L: Length [m]
#        
#        """
#        self.length=length
#        super(oeDrift, self).__init__()
#        
#        self.description=description
#    
#    def propagate(self, wfr, propagation_parameters):
#        """ Overloaded propagate """
#
#        
#        print('Propagating through drift: ' + self.description)
#       
#        beamline = wpg.srwlib.SRWLOptC([], propagation_parameters)
#        srwl.PropagElecField(wfr._srwl_wf, beamline)
           
        
class oeScreen(Empty):
    """
    class: Extends the WPG Screen optical element
    to write image file to common format that is convenient for inspection
    """
    
    def __init__(self, filename=None, 
                       writeIntensity = True,
                       writePhase = False,
                       showIntensity = False,
                       showPhase = False,
                       calculateMetrics = True,
                       description = None,
                       ROI=None):
        """ Constructor for the Screen class.

        :param filename: Name of file to store wavefront data.
        :type filename: str
        :raise IOError: File exists.
        """
        
         # Initialize base class.
        super(oeScreen, self).__init__()
        self.description = description
        
        if filename is not None:
            self.extension = '.h5' #os.path.splitext(filename)[1]
            self.filename = os.path.splitext(filename)[0]
        else:
            self.filename = None
        
        self.writeIntensity = writeIntensity
        self.writePhase = writePhase
        
        self.showIntensity = showIntensity
        self.showPhase = showPhase
        
        self.calculateMetrics = calculateMetrics
        self.ROI = ROI
        
    def describe(self,wfr,ROI):
        ''' 
        print and return summary of wfr
        
        param ROI: tuple of tuples: [{roi_xmin,roi_xmax}, {roi_ymin,roi_ymax}] where each value is in pixels
        '''
              
        #mesh = wfr.params.Mesh           
        [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wfr)
        dx = (xmax-xmin)/nx
        dy = (ymax-ymin)/ny

        intensity = wfr.get_intensity(polarization='horizontal')
        intensityOnAxis = intensity[nx//2, ny//2] / intensity.max()
        
        if ROI is not None:    
            # Note that ROI valus are in units of pixels
            print ('Analysing ROI: [{},{}], [{},{}]'.format(ROI[0][0],ROI[0][1],ROI[1][0],ROI[1][1]))
            intensity = intensity[ROI[0][0]:ROI[0][1],ROI[1][0]:ROI[1][1]]
            xmin, xmax = ROI[0][0]*dx, ROI[0][1]*dx
            ymax, ymin = ROI[1][0]*dy, ROI[1][1]*dy
            nx,ny = np.abs(ROI[0][1]-ROI[0][0]), np.abs(ROI[1][0]-ROI[1][1])
            
        intensityTotal = intensity.sum()
        #intensityTotal = np.squeeze(intensity)
        intensityTotal = intensityTotal*dx*dy#*1e6*1e-9   # [GW] /  dx [mm]/ dy [mm], i.e. GW/mm^2
        
        energy = intensity*wfr.params.photonEnergy/J2EV#*1e3
        try:
            imax = np.max(energy)
        except:
            imax = 0
        
        x_center = intensity[nx//2, :]
        fwhm_x = len(x_center[x_center > x_center.max()/2])*dx

        y_center = intensity[:, ny//2]
        fwhm_y = len(y_center[y_center > y_center.max()/2])*dy
        
        intensityPeak = imax*1e-9*1e6*2*np.pi*(fwhm_x/2.35)*(fwhm_y/2.35)
        
        print('stepX, stepY [um]:', dx * 1e6, dy * 1e6, '\n')
        print('Total power: %g [GW]' % intensityTotal)
        print('Peak power calculated using FWHM:         %g [GW]' %(intensityPeak))
        print('Max irradiance: %g [GW/mm^2]'    %(imax*1e-9))
        
        #label4irradiance = 'Irradiance (W/$mm^2$)'

        summary = {'IntensitySum':intensityTotal,
                   'IntensityPeak':intensityPeak,
                   'dx':dx,
                   'dy':dy,
                   'fwhm_x':fwhm_x, 
                   'fwhm_y': fwhm_y}
        
        print(summary)
        
        return summary

    
    def show(self,wfr):
        
        if self.showIntensity:
            plotWavefront(wfr, self.description)
            #plot_intensity_map(wfr, save=self.description, range_x=None, range_y=None,im_aspect='equal')
            
        if self.showPhase:
            
            plotWavefront(wfr, self.description, phase=True)
            
        
        
    def write(self,wfr):
        
        if self.filename is not None:
            if self.extension == '.h5':
                wfr.store_hdf5(self.filename + self.extension)
                
            if self.writeIntensity == True:
                imageio.imwrite(self.filename + '_intensity.tif', wfr.get_intensity(polarization='horizontal'))
                
            if self.writePhase == True:
                imageio.imwrite(self.filename + '_phase.tif', wfr.get_phase(polarization='horizontal'))
        else:
            print('No filename.  Wavefield at screen will not be saved')
         
        
    def propagate(self, wfr):
        """ Overloaded propagation for this element. """

        #we could propagate to allow resizing etc, but then we need to manage preservation of original wavefield.
        # better to keep functinoality restricted to that of a screen, and manage its use outside class.
        ##super(oeScreen, self).propagate(wfr, propagation_parameters)
        
        if self.filename is not None:
            self.write(wfr)
        
        self.show(wfr)
        
        if self.calculateMetrics:
            metrics = self.describe(wfr, self.ROI)
            return metrics
        else:
            return None
                 
            
    def plotPhase(self, wfr, title, cuts=False, interactive=True):
    #draw wavefront with common functions
        J2EV = 6.24150934e18
        
        wf_phase = wfr.get_phase(polarization='horizontal')
        ii = wf_phase
   
        ii = ii*wfr.params.photonEnergy/J2EV#*1e3
        imax = np.max(ii)
        [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(wfr)
        dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
        print('stepX, stepY [um]:', dx * 1e6, dy * 1e6, '\n')

        if wfr.params.wEFieldUnit != 'arbitrary':
            print('Total power (integrated over full range): %g [GW]' %(ii.sum(axis=0).sum(axis=0)*dx*dy*1e6*1e-9))
            print('Max irradiance: %g [GW/mm^2]'    %(imax*1e-9))
            label4irradiance = 'Irradiance (W/$mm^2$)'
        else:
                ii = ii / imax
                label4irradiance = 'Irradiance (a.u.)'

        [x1, x2, y1, y2] = wfr.get_limits()
    
       
        pylab.figure(figsize=(21,6))
        pylab.imshow(ii, extent=[x1 * 1e3, x2 * 1e3, y1 * 1e3, y2 * 1e3])
        pylab.set_cmap('hot')
        pylab.axis('tight')
        #pylab.colorbar(orientation='horizontal')
        pylab.xlabel('x (mm)')
        pylab.ylabel('y (mm)')
        pylab.axes().set_aspect(0.5)
    
        pylab.title(title)
        pylab.show()

        


class oeResample(WPGOpticalElement):

    """Optical element: Empty.
    This is empty propagator used for sampling and zooming wavefront
    """

    def __init__(self, propagation_parameters, description=None):
        super(oeResample, self).__init__()
        self.propagation_parameters = propagation_parameters
        self.description = description

    def __str__(self):
        return self.description

    def propagate(self, wfr, propagation_parameters):
        """
        Propagate wavefront through empty propagator,
        used for sampling and resizing wavefront
        """
        beamline = wpg.srwlib.SRWLOptC([], self.propagation_parameters)
        srwl.PropagElecField(wfr._srwl_wf, beamline)





def plot_map_I(wf, save='', range_x=None, range_y=None,im_aspect='equal'):
    """
    Plot wavefront in  R-space.

    :param wf: wavefront structure
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
    :param im_aspect: aspect for 2D image, string or float number, see matplotlib set_aspect().
    """
    import matplotlib.pyplot as plt
    import numpy
    # Get the wavefront and integrate over time.
    wf_intensity = wf.get_intensity().sum(axis=-1)

    # Get average and time slicing.
    #average = averaged_intensity(wf, bPlot=True)
    nslices = wf.params.Mesh.nSlices
    if (nslices>1):
        dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)/(nslices-1)
        t0 = dt*nslices/2 + wf.params.Mesh.sliceMin
    else:
        t0 = (wf.params.Mesh.sliceMax+wf.params.Mesh.sliceMin)/2

    # Setup a figure.
    figure = plt.figure(figsize=(10, 10), dpi=100)
    plt.axis('tight')
    # Profile plot.
    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    # Get limits.
    xmin, xmax, ymax, ymin = wf.get_limits()

    # Plot profile as 2D colorcoded map.
    profile.imshow(
        wf_intensity, extent=[xmin*1.e3, xmax*1.e3, ymax*1.e3, ymin*1.e3], cmap="YlGnBu_r")

   # profile.set_aspect(im_aspect, 'datalim')


    # Get x and y ranges.
    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    #x = numpy.linspace(xmin*1.e3,xmax*1.e3,wf_intensity.shape[0])
    #y = numpy.linspace(ymin*1.e3,ymax*1.e3,wf_intensity.shape[1])
    x = numpy.linspace(xmin*1.e3, xmax*1.e3, wf_intensity.shape[1])
    y = numpy.linspace(ymin*1.e3, ymax*1.e3, wf_intensity.shape[0])

    # Labels.
    profile.set_xlabel('$mm$', fontsize=12)
    profile.set_ylabel('$mm$', fontsize=12)

    # x-projection plots above main plot.
    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    print(x.shape, wf_intensity.sum(axis=0).shape)

    x_projection.plot(x, wf_intensity.sum(axis=0), label='x projection')

    # Set range according to input.
    if range_x is None:
        profile.set_xlim([xmin*1.e3, xmax*1.e3])
    else:
        profile.set_xlim([-range_x/2., range_x/2.])

    # Set title.
    if(wf.params.wDomain=='time'):
        x_projection.set_title('t0={:03.1g} s '.format(t0))
    else: #frequency domain
        x_projection.set_title('E0={:05.2g} eV'.format(t0))

    # y-projection plot right of main plot.
    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label='y projection')

    # Hide minor tick labels, they disturb here.
    plt.minorticks_off()

    # Set range according to input.
    if range_y is None:
        profile.set_ylim([ymin*1.e3, ymax*1.e3])
    else:
        profile.set_ylim([-range_y/2., range_y/2.])

    # If requested, save to disk, otherwise show in interactive window.
    if save != '':
        # Add parameters.
        plt.savefig(save)
    else:
        plt.show()
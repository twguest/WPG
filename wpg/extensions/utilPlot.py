#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 22:29:45 2019

@author: gvanriessen
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import matplotlib.ticker as ticker


font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

from math import trunc
#plt.style.use('seaborn-white')
import cplot
import numexpr as ne
from numpy import pi

def update_label(old_label, exponent_text):
    if exponent_text == "":
        return old_label
    
    try:
        units = old_label[old_label.index("[") + 1:old_label.rindex("]")]
    except ValueError:
        units = ""
    label = old_label.replace("[{}]".format(units), "")
    
    exponent_text = exponent_text.replace("\\times", "")
    
    return "{} [{} {}]".format(label, exponent_text, units)
    
def format_label_string_with_exponent(ax, axis='both'):  
    """ Format the label string with the exponent from the ScalarFormatter """
    ax.ticklabel_format(axis=axis, style='sci')

    axes_instances = []
    if axis in ['x', 'both']:
        axes_instances.append(ax.xaxis)
    if axis in ['y', 'both']:
        axes_instances.append(ax.yaxis)
    
    for ax in axes_instances:
        
        ax.major.formatter._useMathText = True
        plt.draw() # Update the text
        exponent_text = ax.get_offset_text().get_text()
        print(exponent_text)
        label = ax.get_label().get_text()
        print(label)
        
        ax.offsetText.set_visible(False)
        ax.set_label_text(update_label(label, exponent_text))
        



def cb(mmin, mmax,pmin,pmax,n):
    
    phase = np.linspace(pmin, pmax,n)
    mag = np.linspace(mmin, mmax,n)
    
    xx = np.array([phase,]*n)
    yy = np.array([mag,]*n)
     
    #xv,yv = np.meshgrid(xx,yy)
    
    z = np.empty([n,n],dtype=complex)
    for i in range(n):
        for j in range(n):
            #z[i,j] = complex(xv[i,j],yv[i,j])
            z[i,j] = yy[i,j]*np.exp(1j*xx[i,j])
   
    return z

def colorbar(mmin=0, mmax=1,pmin=-pi,pmax=+pi,n=100):
    
    assert mmax > mmin
    assert pmax > pmin
    
    abs_scaling=lambda r: r / (r + 1)
    
    val = cb(mmin, mmax,pmin,pmax,n)
    nx, ny = n, n
      
    hx = (mmax - mmin) / nx
    x = np.linspace(mmin + hx / 2, mmax - hx / 2, nx)
    hy = (pmax - pmin) / ny
    y = np.linspace(pmin + hy / 2, pmax - hy / 2, ny)

    angle = np.arctan2(val.imag, val.real)
    absval_scaled = abs_scaling(np.abs(val)) # np.abs(val)/np.max(np.abs(val))#
    srgb_vals = cplot.main.get_srgb(angle, absval_scaled)

    plt.imshow(
        srgb_vals,
        extent=(x.min(), x.max(), y.max(), y.min()),
        interpolation="nearest",
        origin="lower",
        aspect="equal",
    )
    

def cbfmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def tickformat(x):
    if int(x) == float(x):
        return str(int(x))
    else:
        return str(x)        


def tickformat2(value):
    # Truncate and render as int
    return '{:d}'.format(trunc(value))
    

def plotTrThickness(trOpt): 
    
    import matplotlib.pyplot as plt 
    import matplotlib.ticker as mticker
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    import matplotlib.ticker as ticker
    from extensions.utilPlot import cbfmt
    import colorio
    from extensions.utilPlot import cb
    from matplotlib.ticker import FuncFormatter
    
    thickness  = trOpt.getThickness()

    cbLabel = 'Thickness [nm]'
    mesh = trOpt.getSRWTransmissionFunction().mesh
    nx, ny = mesh.nx, mesh.ny
    xmin, xmax = mesh.xStart, mesh.xFin
    ymin, ymax = mesh.yStart, mesh.yFin

    # attempt kludge to use same colormap as for complex value versions of this function
    a = thickness
    b = np.zeros_like(a)
    val = ne.evaluate("complex(a,b)")
    
    hx = (xmax - xmin) / nx
    x = np.linspace(xmin + hx / 2, xmax - hx / 2)
    hy = (ymax - ymin) / ny
    y = np.linspace(ymin + hy / 2, ymax - hy / 2)
     
    angle = np.arctan2(val.imag, val.real)
    absval_scaled = np.abs(val) / (np.max(np.abs(val)))  
    srgb_vals = cplot.main.get_srgb(angle, absval_scaled)

    #setup axes with gridspec
    fig = plt.figure(figsize=(8,6))
    grid = plt.GridSpec(4,4, hspace=0.3, wspace=0.3)
    main_ax = fig.add_subplot(grid[0:2,1:3])
    #y_prof = fig.add_subplot(grid[0:2,0:1],xticklabels=[],sharey=main_ax)
    #x_prof = fig.add_subplot(grid[2,1:3], yticklabels=[], sharex=main_ax)
    cbf = fig.add_subplot(grid[0:2,3])#, yticklabels=[], xticklabels=[])    
    cbf.set_ylabel(cbLabel)
      
    im_main = main_ax.imshow( srgb_vals,
                              extent=(x.min(), x.max(), y.max(), y.min()),
                              interpolation="nearest",
                              origin="lower",
                              aspect="equal",
                             )
    #main_ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True,useOffset=False))
   
    
    #main_ax.ticklabel_format(axis='both',style='sci',useMathText=True)    
    #format_label_string_with_exponent(main_ax, axis='both')
     
    #main_ax.tick_params(axis='both', which='major')
    #main_ax.yaxis.set_major_formatter(EpiCycleScalarFormatter())
    #main_ax.xaxis.set_major_formatter(EpiCycleScalarFormatter())
    #plt.ticklabel_format(style='sci', axis='both', scilimits=(0,3))

    
    main_ax.tick_params(axis='both', which='major')
    fmtx = FuncFormatter(lambda x, pos: tickformat2(x * 1.e6))
    fmty = FuncFormatter(lambda y, pos: tickformat2(y * 1.e6))
    main_ax.xaxis.set_major_formatter(fmtx)
    main_ax.yaxis.set_major_formatter(fmty)

    #main_ax.set_xlabel('x [$10^{-6}$ m]')
    main_ax.set_xlabel('x [$\mu$m]')
    main_ax.set_ylabel('y [$\mu$m]')
    
    mmin, mmax = np.min(thickness), np.max(thickness)    
    nc=100  # number of color levels
    cval = cb(mmin, mmax,0,0,nc)
    
    x = np.linspace(1, 10,1)
    hy = (mmax - mmin) / nc
    y = np.linspace(mmin + hy / 2, mmax - hy / 2, nc)

    angle = np.arctan2(cval.imag, cval.real)
    absval_scaled = np.abs(cval) / np.max(np.abs(cval))
    srgb_cvals = cplot.main.get_srgb(angle, absval_scaled)
    srgb_cvals  =  np.rot90(srgb_cvals) #rotate 90 degrees
        
    cbf.imshow(
         srgb_cvals,
         interpolation="none",
         origin="lower",
         aspect=5,
         extent = [np.min(thickness)*1.e9,np.max(thickness)*1.e9,np.min(thickness)*1.e9,np.max(thickness)*1.e9]
               
    )
    cbf.get_xaxis().set_ticks([])
    
    plt.show()


def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    
    

def plotTrPhaseShift(trOpt,wavelength): 
    

    import matplotlib.pyplot as plt 
    import matplotlib.ticker as mticker
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from extensions.utilPlot import cbfmt
    import colorio
    from utilPlot import cb
    
    phase = trOpt.getPhaseShift(wavelength)

    cbLabel = 'Phase shift [rad.]'
    mesh = trOpt.getSRWTransmissionFunction().mesh
    nx, ny = mesh.nx,  mesh.ny
    xmin, xmax = mesh.xStart, mesh.xFin
    ymin, ymax = mesh.yStart, mesh.yFin

    # attempt kludge to use same colormap as for complex value versions of this function
    
    b = phase
    a = np.zeros_like(b)
    val = ne.evaluate("complex(a,b)")
    
    hx = (xmax - xmin) / nx
    x = np.linspace(xmin + hx / 2, xmax - hx / 2)
    hy = (ymax - ymin) / ny
    y = np.linspace(ymin + hy / 2, ymax - hy / 2)

    angle = np.arctan2(val.imag, val.real)
    absval_scaled = np.abs(val) / (np.max(np.abs(val)))
    srgb_vals = cplot.main.get_srgb(angle, absval_scaled)

    #setup axes with gridspec
    fig = plt.figure(figsize=(8,6))
    grid = plt.GridSpec(4,4, hspace=0.3, wspace=0.3)
    main_ax = fig.add_subplot(grid[0:2,1:3])
    #y_prof = fig.add_subplot(grid[0:2,0:1],xticklabels=[],sharey=main_ax)
    #x_prof = fig.add_subplot(grid[2,1:3], yticklabels=[], sharex=main_ax)
    cbf = fig.add_subplot(grid[0:2,3])#, yticklabels=[], xticklabels=[])    
    cbf.set_ylabel(cbLabel)
      
    im_main = main_ax.imshow( srgb_vals,
                              extent=(x.min(), x.max(), y.max(), y.min()),
                              interpolation="nearest",
                              origin="lower",
                              aspect="equal",
                             )
    #main_ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True,useOffset=False))
   
    
    #main_ax.ticklabel_format(axis='both',style='sci',useMathText=True)    
    #format_label_string_with_exponent(main_ax, axis='both')
     
    

    #main_ax.tick_params(axis='both', which='major')
    #main_ax.yaxis.set_major_formatter(EpiCycleScalarFormatter())
    #main_ax.xaxis.set_major_formatter(EpiCycleScalarFormatter())
    #plt.ticklabel_format(style='sci', axis='both', scilimits=(0,3))
    
    
    main_ax.tick_params(axis='both', which='major')
    fmtx = FuncFormatter(lambda x, pos: tickformat2(x * 1e6))
    fmty = FuncFormatter(lambda y, pos: tickformat2(y * 1e6))
    main_ax.xaxis.set_major_formatter(fmtx)
    main_ax.yaxis.set_major_formatter(fmty)

    #main_ax.set_xlabel('x [$10^{-6}$ m]')
    main_ax.set_xlabel('x [$\mu$m]')
    main_ax.set_ylabel('y [[$\mu$m]')
    
    mmin, mmax = np.min(thickness), np.max(thickness)    
    nc=100  # number of color levels
    cval = cb(mmin, mmax,0,0,nc)
    
    x = np.linspace(1, 10,1)
    hy = (mmax - mmin) / nc
    y = np.linspace(mmin + hy / 2, mmax - hy / 2, nc)

    angle = np.arctan2(cval.imag, cval.real)
    absval_scaled = np.abs(cval) / np.max(np.abs(cval))
    srgb_cvals = cplot.main.get_srgb(angle, absval_scaled)
    #srgb_cvals  =  np.rot90(srgb_cvals) #rotate 90 degrees
        
    cbf.imshow(
         srgb_cvals,
         interpolation="none",
         origin="lower",
         aspect=5,
        # extent = [np.min(thickness)*1.e9,np.max(thickness)*1.e9,np.min(thickness)*1.e9,np.max(thickness)*1.e9]
               
    )
    cbf.get_xaxis().set_ticks([])
    
    plt.show()


def plotWavefieldComplex(wf,  abs_scaling=lambda r: r / (r + 1)):
    
    mesh = wf.params.Mesh
    nx, ny = mesh.nx,  mesh.ny
    xmin, xmax = mesh.xMin, mesh.xMax
    ymin, ymax = mesh.yMin, mesh.yMax

    b = wf.get_imag_part(slice_number=0)
    a = wf.get_real_part(slice_number=0)
    val = ne.evaluate("complex(a,b)")
    
    
    hx = (xmax - xmin) / nx
    x = np.linspace(xmin + hx / 2, xmax - hx / 2, nx)
    hy = (ymax - ymin) / ny
    y = np.linspace(ymin + hy / 2, ymax - hy / 2, ny)

    angle = np.arctan2(val.imag, val.real)
    angle = (angle-np.min(angle))   / np.max(angle-np.min(angle))  
    
    mag =np.abs(val)
    mag = (mag-np.min(mag))/np.max(mag-np.min(mag))
    #mag = abs_scaling(mag)                                                                                                                                                                                        
    
    srgb_vals = cplot.main.get_srgb(angle, mag)
    
    #setup axes with gridspec
    fig = plt.figure(figsize=(8,6))
    grid = plt.GridSpec(4,4, hspace=0.3, wspace=0.3)
    main_ax = fig.add_subplot(grid[0:2,1:3])
    #y_prof = fig.add_subplot(grid[0:2,0:1],xticklabels=[],sharey=main_ax)
    #x_prof = fig.add_subplot(grid[2,1:3], yticklabels=[], sharex=main_ax)
    cbf = fig.add_subplot(grid[0:2,3])#, yticklabels=[], xticklabels=[])
    cbf.set_xlabel('mag.')
    cbf.set_ylabel('phase (rad.)')
    
    im_main = main_ax.imshow( srgb_vals,
                              extent=(x.min(), x.max(), y.max(), y.min()),
                              interpolation="none",
                              origin="lower",
                              aspect="equal"
                             )
    
    main_ax.tick_params(axis='both', which='major')
    #main_ax.set_xlabel('x [$10^{-6}$ m]')
    main_ax.set_xlabel('x [$\mu$m]')
    main_ax.set_ylabel('y [[$\mu$m]')
    
    fmtx = FuncFormatter(lambda x, pos: tickformat2(x * 1.e6))
    fmty = FuncFormatter(lambda y, pos: tickformat2(y * 1.e6))
    main_ax.xaxis.set_major_formatter(fmtx)
    main_ax.yaxis.set_major_formatter(fmty)

    
    # prepare and format colorbar
    #min=0, mmax=1, pmin=-math.pi, pmax=+math.pi,
    mmin,mmax=np.min(mag),np.max(mag)
    pmin,pmax=np.min(angle),np.max(angle)
    #pmin=-pi
    #pmax=+pi
    
    print('Phase: [%3.3f,%3.3f]' % (pmin,pmax))
    print('Abs: [%3.3f,%3.3f]' % (mmin,mmax))
    
    n=100
    cval = cb(mmin, mmax,pmin,pmax,n)                                                                          
    nx, ny = n,n
      
    hx = (mmax - mmin) / nx
    x = np.linspace(mmin + hx / 2, mmax - hx / 2, nx)
    hy = (pmax - pmin) / ny
    y = np.linspace(pmin + hy / 2, pmax - hy / 2, ny)

    angle = np.arctan2(cval.imag, cval.real)
    abscval_scaled = abs_scaling(np.abs(cval))
    srgb_cvals = cplot.main.get_srgb(angle, abscval_scaled)
    
    cbf.imshow(
        srgb_cvals,
        extent=(x.min(), x.max(), y.max(), y.min()),
        interpolation="nearest",
        origin="lower",
        aspect=0.25
    )
    #cbf.grid(color='w', linestyle='-',linewidth=1)
    
    cbf.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 4))
    cbf.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 8))
    cbf.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

    cbf.xaxis.set_ticks([mmin, (mmax-mmin) / 2., mmax])
    forceAspect(cbf,0.4)

    # plt.show()
    
    return val



def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
#https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib    
    
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\dfrac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\dfrac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\dfrac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
        
    def locator(self):
            return plt.MultipleLocator(self.number / self.denominator)
        
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
    
    

def valid_imshow_data(data):
    data = np.asarray(data)
    if data.ndim == 2:
        return True
    elif data.ndim == 3:
        if 3 <= data.shape[2] <= 4:
            return True
        else:
            print('The "data" has 3 dimensions but the last dimension '
                  'must have a length of 3 (RGB) or 4 (RGBA), not "{}".'
                  ''.format(data.shape[2]))
            return False
    else:
        print('To visualize an image the data must be 2 dimensional or '
              '3 dimensional, not "{}".'
              ''.format(data.ndim))
        return False
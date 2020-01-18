#!/usr/bin/env python

"""
This script is used for plotting complex data

I created some, but some functions were also stolen from Ptypy
Mine generally produce nicer plots for some reason

"""

import os
import sys
import h5py
import glob
from glob import iglob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar.scalebar import SI_LENGTH
from scipy import ndimage
from matplotlib import ticker
import matplotlib.font_manager as fm
import matplotlib

import matplotlib.patches as patches

from matplotlib.colors import hsv_to_rgb









def Complex_colour_old(complex_array, pixelSizeAtSample, title, savepath, type):
    from matplotlib.colors import hsv_to_rgb
    import numpy as np
    import pylab as plt
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib_scalebar.scalebar import SI_LENGTH
    
    
    # type changes positioning of colour scale plot thing
    def normalize(M):
        return (M-np.min(M))/(np.max(M)-np.min(M))    
    
    phase = np.angle(complex_array)
    magnitude = np.abs(complex_array)
    
    phase = np.swapaxes(phase,0,1)
    magnitude = np.swapaxes(magnitude,0,1)
    
    # Pre-defining the rgb array size
    rgbArray = np.zeros((int(np.shape(magnitude)[0]),int(np.shape(magnitude)[1]),3), 'uint8')
    
    # Normalising the phase and magnitude values
    Norm_phase = normalize(phase)*2*np.pi
    Norm_magnitude = normalize(magnitude)
    
    # Defining the rgb channels from the normalised phase and magnitude values
    r = 0.5*(np.sin(Norm_phase)+1)*Norm_magnitude*255
    g = 0.5*(np.sin(Norm_phase+np.pi/2)+1)*Norm_magnitude*255
    b = 0.5*(-np.sin(Norm_phase)+1)*Norm_magnitude*255
            
    # making the rgb array of the complex field
    rgbArray[...,0] = r
    rgbArray[...,1] = g
    rgbArray[...,2] = b
            
    # Creating the colourmap
    V, H = np.mgrid[0:1:100j, 0:1:300j]
    S = 0.7*np.ones_like(V)
    HSV = np.dstack((H,S,V))
    RGB = hsv_to_rgb(HSV)
    
    color = RGB.reshape((RGB.shape[0]*RGB.shape[1],RGB.shape[2]))
    theta, R = np.meshgrid(
        np.linspace(0,2*np.pi,301),
        np.linspace(0,1,100),
    )
    
    t,r = np.meshgrid(
        np.linspace(0,1,220),
        np.linspace(0,1,130),
    )    
    

    # Plotting the rgbArray
    fig,ax = plt.subplots()
    ax.imshow(rgbArray)
    ax.set_position([0.0, 0.0, 1, 1])
    ax.axis('off')
    #ax.set_title(
    #        title, color='black'#, fontweight='#from printline import printlinrbold'
    #)
    scalebar = ScaleBar(pixelSizeAtSample, 'm', SI_LENGTH, box_alpha=1)
    ax.add_artist(scalebar)
    ax.text(0.01,0.85, title, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='white', fontsize=18, fontweight='bold')
    
    # Adding the colourmap
    if type==0:
        ax2 = fig.add_axes([0.75, 0.1, 0.2, 0.2],projection='polar', facecolor='black') #plotting the colourmap in the bottom right corner
    elif type ==1:
        ax2 = fig.add_axes([0.6, 0.1, 0.2, 0.2],projection='polar', facecolor='black') #plotting the colourmap in the bottom right corner
    ax2.pcolormesh(
        theta,R,
        np.zeros_like(R),
        color = color,
    )
    
    ax2.set_thetagrids(np.linspace(0,2*np.pi,5)[:-1], frac=1.15)
    ax2.set_xticks(np.linspace(0,2*np.pi,5)[:-1])
    #ax2.set_xticks([ np.pi, 2*np.pi,0][:-1])
    ax2.set_xticklabels(
       #['%0.2f        ' % (np.min(phase)),'       %0.2f' % (np.max(phase))] , color='white', fontweight='bold'
       #[r'$\boldsymbol{\pi}$        ' ,r'       $\boldsymbol{0/2\pi}$'] , color='white', fontweight='bold'
       ['',r'$2\pi$ / $0$', '' ,r'$\pi$',''] , color='white', fontweight='bold', fontsize=12
    )
    
    ax2.set_yticks(np.linspace(0,1,2))
    
    ax2.set_yticklabels(
        #['{:.2}'.format(i) for i in np.linspace(np.min(magnitude),np.max(magnitude),2)], color='white', fontweight='bold'
        ['0','1'], color='white', fontsize=12
    )
    ax2.grid('off')
    # shamelessly taken from stack overflow:
    tick = [ax2.get_rmax(),ax2.get_rmax()*0.95]
    for t  in np.deg2rad(np.arange(0,360,45)):
        ax2.plot([t,t], tick, lw=0.72, color="k")
    ax2.plot((0,0.4), (0,1), c='k', linewidth=0.8)
    plt.savefig(savepath)
    plt.close("all") # This needs to be repeated throughout script to reduce memory. Computer will get angry otherwise






def Colour_complex(complex_array,
				   ImageOut=0,
                   norm_phase=False,
                   norm_mag=False,
                   smooth_complex=False,
                   add_colourmap=False,
                   add_scalebar= [False, 0],
                   showfigs=0,
                   clim_mag = 0,
                   clim_phase = 0,
                   save_solo = 0):


	
	
	# Extracting the magnitude and phase components	
	phase = np.angle(complex_array)
	magnitude = np.sqrt(abs(complex_array))
	
	# Smoothing the magnitude and phase. Need to do individually
	if smooth_complex==True:
		phase = Smooth_array(phase)
		magnitude = Smooth_array(magnitude)
	
	
	# Fitting the image to fit within a certain range. Possibly used to 
	# be consistent with other images. Easiest way is to alter corners of
	# image to be the max and min values in clim
	if type(clim_mag) is tuple:
		magnitude[0,0] = clim_mag[0]
		magnitude[-1,-1] = clim_mag[1]
	
	if type(clim_phase) is tuple:
		phase[0,0] = clim_phase[0]
		phase[-1,-1] = clim_phase[1]
	
	
		
	
	
	# Normalising the phase and magnitude values
	# This used to stretch values out over whole range, enhances contrast in final image
	if norm_phase==True:
		phase = normalize(phase)*2*np.pi
		
	if norm_mag==True:
		magnitude = normalize(magnitude)

	# Pre-defining the rgb array size
	rgbArray = np.zeros((int(np.shape(magnitude)[0]),int(np.shape(magnitude)[1]),3), 'uint8')
	
	
	# Defining the rgb channels from the normalised phase and magnitude values
	r = 0.5*(np.sin(phase)+1)*magnitude*255
	g = 0.5*(np.sin(phase+np.pi/2)+1)*magnitude*255
	b = 0.5*(-np.sin(phase)+1)*magnitude*255

	# making the rgb array of the complex field
	rgbArray[...,0] = r
	rgbArray[...,1] = g
	rgbArray[...,2] = b
 

	# Plotting the rgbArray
	fig,ax = plt.subplots()
	ax.imshow(rgbArray, origin='lower')
	#ax.set_position([0.0, 0.0, 1, 1])
	ax.set_position([0.1, 0.1, 0.7, 0.7])
	ax.axis('off')
	
	if add_scalebar[0] == True:
		scalebar = ScaleBar(add_scalebar[1], 'm', SI_LENGTH, box_alpha=0, location=9, length_fraction=0.6, height_fraction = 0.06, frameon='True', border_pad=-3, pad=0)
		plt.gca().add_artist(scalebar)
	
	if add_colourmap==True:
		
		# Creating the colourmap
		V, H = np.mgrid[0:1:100j, 0:1:300j]
		S = 1*np.ones_like(V)
		HSV = np.dstack((H,S,V))
		RGB = hsv_to_rgb(HSV)
		color = RGB.reshape((RGB.shape[0]*RGB.shape[1],RGB.shape[2]))


		theta, R = np.meshgrid(
			np.linspace(0,2*np.pi,301),
			np.linspace(0,1,100),
		)


		t,r = np.meshgrid(
			np.linspace(0,1,220),
			np.linspace(0,1,130),
		)
		
		# Adding the colourmap
		ax2 = fig.add_axes([0.8, 0.2, 0.2, 0.2],projection='polar', facecolor='black') #plotting the colourmap in the bottom right corner
		ax2.pcolormesh(
			theta,R,
			np.zeros_like(R),
			color = color,
		)

		ax2.set_thetagrids(np.linspace(0,2*np.pi,5)[:-1], frac=1.15)
		ax2.set_xticks(np.linspace(0,2*np.pi,5)[:-1])

		if norm_phase==True:
			ax2.set_xticklabels(['',r'$2\pi$ / $0$', '' ,r'$\pi$',''] , color='black', fontweight='bold', fontsize=12)
		else:
			low_bound = np.min(phase) #+ (0.001*abs(np.min(phase)))
			mid_bound = (np.max(phase)+np.min(phase))/2
			high_bound = np.max(phase) #- (0.001*abs(np.max(phase)))
			ax2.set_xticklabels(['',r'$'+str(high_bound)+'$ / $'+str(low_bound)+'$', '' ,r'$'+str(mid_bound)+'$',''] , color='white', fontweight='bold', fontsize=12, format='%0.2f')
		ax2.set_yticks(np.linspace(0,1,2))
		if norm_mag==True:
			ax2.set_yticklabels(['0','1'], color='white', fontsize=12)
		else:
			ax2.set_yticklabels([str(np.min(magnitude)),str(np.min(magnitude))], color='white', fontsize=12, format='%0.2f')
		ax2.grid('off')
		
		ticklabels = [t for t in plt.gca().get_yticklabels()]
		ticklabels[1].set_color("black")

		tick = [ax2.get_rmax(),ax2.get_rmax()*0.95]
		for t  in np.deg2rad(np.arange(0,360,45)):
			ax2.plot([t,t], tick, lw=0.72, color="k")
		ax2.plot((0,0.4), (0,1), c='k', linewidth=0.8)

	
	
	if type(ImageOut)==str:
		print('saving file to: ' + str(ImageOut) + '_coloured.png')
		plt.savefig(ImageOut + '_coloured.png', transparent=True)
	
	if showfigs:
		plt.show()
	
	
	ax1 = plt.gca()
	pos_ax1 = ax1.get_position()
	if save_solo:
		plt.figure(20)
		plt.imshow(rgbArray, extent=None, interpolation='none', origin='lower')
		plt.axis('off')
		
		
		ax2 = plt.gca()
		ax2.set_position(pos_ax1)
		if showfigs:
			plt.show()
		plt.savefig(ImageOut +'_colored_solo.png', transparent=True, bbox_inches='tight', pad_inches=0)
	
	
	plt.close('all')
	
	return rgbArray


	


def complex2hsv(cin, vmin=0., vmax=None, norm_hue=False, norm_value=False, smoothing=False):
	"""\
	Transforms a complex array into an RGB image,
	mapping phase to hue, amplitude to value and
	keeping maximum saturation.

	Parameters
	----------
	cin : ndarray
		Complex input. Must be two-dimensional.

	vmin,vmax : float
		Clip amplitude of input into this interval.

	Returns
	-------
	rgb : ndarray
		Three dimensional output.

	See also
	--------
	complex2rgb
	hsv2rgb
	hsv2complex
	"""
	# HSV channels
	phase = np.angle(cin)
	mag = abs(cin)
	
	
	if smoothing==True:
		h = .5*Smooth_array(phase)/np.pi + .5
		s = np.ones(cin.shape)
		v = Smooth_array(mag)
	else:
		h = .5*phase/np.pi + .5
		s = np.ones(cin.shape)
		v = mag
	

	if norm_hue:
		h = normalize(h)
	
	if norm_value:
		v = normalize(v)
	
	#v = abs(cin)/np.max(abs(cin))
	if vmin is None:
		vmin = v.min()
	if vmax is None:
		vmax = v.max()
	if vmin==vmax:
		v = np.ones_like(v) * v.mean()
		v = v.clip(0.0, 1.0)
	else:
		assert vmin < vmax
		v = (v.clip(vmin, vmax)-vmin)/(vmax-vmin)
	
	return np.asarray((h, s, v))



# I didn't create this one
def complex2rgb(cin, **kwargs):
	"""
	Executes `complex2hsv` and then `hsv2rgb`

	See also
	--------
	complex2hsv
	hsv2rgb
	rgb2complex
	"""
	return hsv2rgb(complex2hsv(cin, **kwargs))



# I didn't create this one
def hsv2rgb(hsv):
	"""\
	HSV (Hue,Saturation,Value) to RGB (Red,Green,Blue) transformation.

	Parameters
	----------
	hsv : array-like
		Input must be two-dimensional. **First** axis is interpreted
		as hue,saturation,value channels.

	Returns
	-------
	rgb : ndarray
		Three dimensional output. **Last** axis is interpreted as
		red, green, blue channels.

	See also
	--------
	complex2rgb
	complex2hsv
	rgb2hsv
	"""
	# HSV channels
	h, s, v = hsv

	i = (6.*h).astype(int)
	f = (6.*h) - i
	p = v*(1. - s)
	q = v*(1. - s*f)
	t = v*(1. - s*(1.-f))
	i0 = (i % 6 == 0)
	i1 = (i == 1)
	i2 = (i == 2)
	i3 = (i == 3)
	i4 = (i == 4)
	i5 = (i == 5)

	rgb = np.zeros(h.shape + (3,), dtype=h.dtype)
	rgb[:, :, 0] = 255*(i0*v + i1*q + i2*p + i3*p + i4*t + i5*v)
	rgb[:, :, 1] = 255*(i0*t + i1*v + i2*v + i3*q + i4*p + i5*p)
	rgb[:, :, 2] = 255*(i0*p + i1*p + i2*t + i3*v + i4*v + i5*q)

	return rgb

def rgb2hsv(rgb):
    """
    Reverse to :any:`hsv2rgb`
    """
    eps = 1e-6
    rgb = np.asarray(rgb).astype(float)
    maxc = rgb.max(axis=-1)
    minc = rgb.min(axis=-1)
    v = maxc
    s = (maxc-minc) / (maxc+eps)
    s[maxc <= eps] = 0.0
    rc = (maxc-rgb[:, :, 0]) / (maxc-minc+eps)
    gc = (maxc-rgb[:, :, 1]) / (maxc-minc+eps)
    bc = (maxc-rgb[:, :, 2]) / (maxc-minc+eps)

    h = 4.0+gc-rc
    maxgreen = (rgb[:, :, 1] == maxc)
    h[maxgreen] = 2.0+rc[maxgreen]-bc[maxgreen]
    maxred = (rgb[:, :, 0] == maxc)
    h[maxred] = bc[maxred]-gc[maxred]
    h[minc == maxc] = 0.0
    h = (h/6.0) % 1.0

    return np.asarray((h, s, v))


def hsv2complex(cin):
    """
    Reverse to :any:`complex2hsv`
    """
    h, s, v = cin
    return v * np.exp(np.pi*2j*(h-.5)) / v.max()


def rgb2complex(rgb):
    """
    Reverse to :any:`complex2rgb`
    """
    return hsv2complex(rgb2hsv(rgb))

if __name__ == "__main__":
    
    from extensions.twpg_wavefront import Wavefront
    from extensions.twpg_uti_wf import fixPhase
    
    import numpy as np
    from numpy import pi
    import pylab as plt
    from colorsys import hls_to_rgb
    from sklearn.preprocessing import minmax_scale as norm
    import seaborn as sns
    
    def colorize(z):
        r = np.abs(z)
        arg = np.angle(z) 
    
        h = (arg + pi)  / (2 * pi) + 0.5
        l = norm(r)
        s = 100
    
        c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
        c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
        c = c.swapaxes(0,2) 
        return c
    



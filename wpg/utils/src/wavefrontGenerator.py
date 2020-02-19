#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/9/19
Upated 18/9/2019.  

@author: Grant van Riessen (La Trobe University)
"""


from wpg.beamline import Beamline
from wpg.wavefront import Wavefront


from wpg.srwlib import SRWLOptD,SRWLOptC,SRWLOptA,SRWLOptT
from extensions.gvrutils import plotWavefront
from extensions.gvrutils import writeIntensity, pixelScale

from wpg.optical_elements import Use_PP


use_gpu = False 

if use_gpu:
    import afnumpy as np
    import afnumpy.fft as fft
    use = 'afnumpy/GPU'
else:
    import numpy as np
    import numpy.fft as fft
use = 'numpy/CPU'

import os

try:
    from wpg import srwlpy as srwl
except ImportError:
    import srwlpy as srwl  #  Hack for read the docs

from wpg.uti_io import *



def gsmSource():
    
    def __init__(self, nModes_x, nModes_y, rho_I, rho_mu, EKeV, dimensions)
        '''
        param: nx - number of horizontal modes
        param: ny - number of vertical modes
        param: EKeV  energy in keV
        param dimensions - list of dimensions [xmin, xmax, ymin, ymax] in units of m
        
        p_mu and p_I are the rms widths of the degree of coherence and of the 
        intensity of the source.
        
        definitions of parameters follow Starikov and Wolf, 1982:

        beta is a measure of the "degree of global coherence" of the source:
                             
                               beta =  p_mu/p_I 
                
        
        When beta >> 1, the source is effectively spatially coherent in the
        global sense and is then found to be well represented by a single mode. 
        When beta << 1, the source is effectively spatially incoherent in the
        global sense, and the number of modes needed to describe its behavior 
        is of the order of beta/3
        
        
        '''
        self.nx = nModes_x
        self.ny = nModes_y
        self.n = self.nx*self.ny
        self.rho_mu = rho_mu
        self.rho_I = rho_I
        self.beta = self.rho_mu / self.rho_I
        
        self.modes2D()
        self.getEigenValues()
        
        self.EkeV = EkeV
        
        
        self.xmin, self.xmax, self.ymin, self.ymax = dimensions
        
    
    def modes2D(self):
        """
        return transverse Gauss-Hermite Mode Order pairs up to order mx = nx and order my=my
    
        """
        rx=range(self.nx)
        ry=range(self.ny)
        self.modes = [(x,y) for x in rx for y in ry]
        self.modes.sort(key=lambda tup: max(tup))
 
    
    def getEigenValues():
        self.eigenValues = [self.eigenVal(i) for i in range(self.n)]
    
    
    def eigenVal(n):
        
        adapt to x, x vals
        
        """
        
        return eigenvalue normalised to eigenvalue of fundamental
        
        definitions follow Starikov and Wolf, 1982:
           
        rho_mu and rho_I are the rms widths of the degree of coherence and of the intensity of the source.
        
        beta =  rho_mu/rho_I is a measure of the "degree of global coherence" of the source.
    
        When beta >> 1, the source is effectively spatially coherent in the global sense and is then
        found to be well represented by a single mode. 
        When beta << 1, the source is effectively spatially incoherent in the global
        sense, and the number of modes needed to describe its behavior is of the order of beta/3
         
        """
        
        a = 1/(2*self.rho_I_x**2)
        b = 1/(2*self.rho_mu_x**2)
        c = (a**2 + 2*a*b)**(1/2)
        l_0=1.0
        
        l_n = l_0*(b/(a+b+c))**range(self.n)
        
        return l_n
    

    def plotEigenValues(e,_threshold=None):
    
        plt.plot(e,'bo')
        if threshold:
            e = self.eigenValues[eigenValues>_threshold]
        else:
            e = self.eigenValues
            
        plt.plot(e,'r+')
        plt.title('Eigenvalues')
        plt.xlabel('Mode')
        plt.ylabel('$\lambda/\lambda_0$')
        plt.legend(loc=2)
        plt.show()
        
        
    def generateMode(self,n, outPath = None):
        ''' 
        Generate the nth mode        
        
        This provides a wrapper for SRW build_gaussian_wavefront_xy
        
        '''
        
        %%%%def generate_mode(self, mode_no, eigenfn, outPath = None):
        
        gauss_0 = primary_gaussian() # properties of gaussian

        gauss_0.d2waist = 0
        print("Generating Mode {}".format(mode_no))
        mx = self.modeIndices[n][0]
        my = self.modeIndices[n][1]
        
        xw = self.rho_I *
        yw = self.rho_I * 
        
        M2x = (2*mx)+1
        M2y = (2*my)+1
        
        wfr = Wavefront(build_gauss_wavefront_xy(gauss_0.nx ,gauss_0.ny,
                                                      self.EkeV,
                                                      self.xMin, self.xMax,
                                                      self.yMin, self.yMax,
                                                      gauss_0.sigX*(M2x**-0.5),
                                                      gauss_0.sigY*(M2y**-0.5),
                                                      gauss_0.d2waist,
                                                      tiltX = 0,
                                                      tiltY = 0,
                                                      _mx = mx,
                                                      _my = my))
        
        if outPath is not None:
            wfr.store_hdf5(outdir + "wfr_mode_{}.hdf5".format(mode_no))
        
        return wfr
        
        
        
        
        
        
        
        
        

def constructTestWaveField(npoints=256,ekeV=0.6,z1=5, show=False):  
    # Simple function for generating  simple Gaussian wavefield for testing purposes


    from wpg.srwlib import srwl_uti_ph_en_conv
    from extensions.gvrutils import calculate_theta_fwhm_cdr
    from wpg.generators import build_gauss_wavefront_xy


    # set wavefield parameters
    qnC        = 0.01
    wlambda    = srwl_uti_ph_en_conv(ekeV, _in_u='keV', _out_u='nm')
    theta_fwhm = calculate_theta_fwhm_cdr(ekeV,qnC)
    k          = 2*np.sqrt(2*np.log(2))
    range_xy   = (theta_fwhm/k*z1*5.)*2.0
    sigX       = 12.4e-10*k/(ekeV*4*np.pi*theta_fwhm)

    # construct wavefield
    wf0=build_gauss_wavefront_xy(nx=npoints, ny=npoints, ekev=ekeV,
                                              xMin=-range_xy/2 ,xMax=range_xy/2,
                                              yMin=-range_xy/2, yMax=range_xy/2,
                                              sigX=sigX, sigY=sigX,
                                              d2waist=z1,
                                              _mx=0, _my=0 )
    wfr = Wavefront(srwl_wavefront=wf0)    
    print('Wavelength=%f, theta FWWM=%f, range XY = %f, sig X = %f' % (wlambda, theta_fwhm, range_xy, sigX))

    wfr._srwl_wf.unitElFld = 1#'sqrt(Phot/s/0.1%bw/mm^2)'

    if show==True: #display the wavefield
        plotWavefront(wfr, 'Wavefield at source')

    return wfr



    
    def generate_mode_test(self, mode_no, eigenfn, outdir, gwaist = None):
        gauss_0 = primary_gaussian() # properties of gaussian
        gauss_0.d2waist = -9.0000
        if gwaist is not None:
            gauss_0.d2waist = gwaist
        print("Generating Mode {}".format(mode_no))
        _mx = eigenfn[mode_no][0]
        _my = eigenfn[mode_no][1]
        
        M2x = (2*_mx)+1
        M2y = (2*_my)+1
        
        GsnBm = SRWLGsnBm(_x = 0, _y = 0, _z = 0,
                          _sigX = 100e-06,#gauss_0.sigX*(M2x**-0.25),
                          _sigY = 100e-06,#gauss_0.sigY*(M2y**-0.25),
                          _xp = 0, #gauss_0.xp,
                          _yp = 0,#gauss_0.yp,
                          _mx = _mx,
                          _my = _my,
                          _avgPhotEn = gauss_0.E)
        
        wfr = Wavefront(build_gaussian(GsnBm,
                                       gauss_0.nx*2, gauss_0.ny*2,
                                       gauss_0.xMin*2, gauss_0.xMax*2,
                                       gauss_0.yMin*2, gauss_0.yMax*2))

        wfr.store_hdf5(outdir + "wfr_mode_{}.hdf5".format(mode_no))
        return wfr
    
        

def test():
    '''
    author: twg
    
    '''
    from joblib import Parallel, delayed, dump, load
    from tqdm import tqdm
    
    nmodes = 30
    eigenval = []
    eigenfn = gen_modes(nmodes, 1)
    
    beta = 0.04258226873188785
    for itr in range(nmodes):
        eigenval.append(calc_eigenval(1, beta, itr))
    eigenval = np.asarray(eigenval)
    eigenfn  = np.asarray(eigenfn)
    np.savetxt(outdir + "eigenfn.csv", eigenfn, delimiter = ",")
    np.savetxt(outdir + "eigenval.csv", eigenval, delimiter = ",")
    
    outdir = r"/nfs/data/users/twg/gsmModes/"
    
    beamline = sxri_bl()
    beamline.setup_OE()
    #beamline.setup_beamline()
    #beamline.setup_metrology_beamline(beamline.bl)
    Parallel(n_jobs = 5,backend= "threading")(delayed(beamline.generate_mode)(mode_no, eigenfn, outdir) for mode_no in tqdm(range(nmodes)))



 def load_mode(self, mode_no, indir):
        wfr = Wavefront()
        wfr.load_hdf5(indir + "/wfr_mode_{}.hdf5".format(mode_no))
        return wfr
    def run_genProp(self, mode_no, eigenfn, outdir):
        print("Initialising Run {}".format(mode_no))
        wfr = self.generate_mode(mode_no, eigenfn)
        self.propagate_metrology_beamline(wfr, mode_no, outdir)
        ###mode.store_hdf5(outdir + "wfr_0.hdf5") 
    
    def run_loadProp(self, mode_no, indir, outdir, atSource = True):
        print("Loading Wavefront: {}".format(mode_no))
        wfr = self.load_mode(mode_no, indir)
        self.propagate_metrology_beamline(wfr, mode_no, outdir, atSource)
        #self.propagate_from_exitslits(wfr, mode_no, outdir)

        
      

    def test_gsm(self, beta = 1, nx = 10, ny = 1, pos = "end", export = False, outdir = ""):
        """
        Initialise GSM beam as a container of multiple coherent modes

        Modes dictionary is setup by order {wfr, eigenvalue/weight}

        :param p_mu: rms width of the degree of coherence
        :param p_I: rms width of the intensity profile
        :param n: number of values for eigenvalue calculation (almost arb.)
        """
        print("Initialising Gaussian-Schell Model Beam")
        self.log.info("Initialising Gaussian-Schell Model Beam")

        gauss_0 = primary_gaussian() # properties of gaussian
        N = nx*ny
        

        results = []

        if export == True:
            os.mkdir(outdir)
        else:
            pass

        for itr in range(N):
            self.wfr.eigenval.append(eigenvaluePartCoh(1, beta,itr))

        self.wfr.eigenval = self.wfr.eigenval[0:N]
        self.wfr.eigenfn = gen_modes(nx, ny)

        if pos == "start":
            gauss_0.d2waist = 0
        elif pos == "end":
            gauss_0.d2waist = 2.012
        else:
            assert("wavefront posn' should be 'start' or 'end'")

        print("Generating {} Coherent Modes".format(N))
        self.log.info("Generating {} Coherent Modes".format(N))
        
        eField = np.zeros((2048, 2048, 1))

        if pos == "start":
            gauss_0.d2waist = 0
        elif pos == "end":
            gauss_0.d2waist = 2.012
        else:
            assert("wavefront posn' should be 'start' or 'end'")
        
        print(self.wfr.eigenfn)
        print("Generating {} Coherent Modes".format(N))
        self.log.info("Generating {} Coherent Modes".format(N))
        for itr in range(N):
            
            res = []
            res.append(self.wfr.eigenfn[itr])
            res.append(self.wfr.eigenval[itr])

            if itr % 5 == 0:
                print("Generating Mode {}/{}".format(itr+1, N))
                self.log.info("Generating Mode {}/{}".format(itr+1, N))
                
            _mx = self.wfr.eigenfn[itr][0]
            _my = self.wfr.eigenfn[itr][1]

            M2x = (2*_mx)+1
            M2y = (2*_my)+1
            self.wfr.modes.append(Wavefront(build_gauss_wavefront_xy(gauss_0.nx ,gauss_0.ny,
                                           self.global_E*1e-03,
                                           gauss_0.xMin, gauss_0.xMax,
                                           gauss_0.yMin,gauss_0.yMax,
                                           gauss_0.sigX*M2x**-0.25,
                                           gauss_0.sigY*M2y**-0.25,
                                           gauss_0.d2waist,
                                           tiltX = 0,
                                           tiltY = 0,
                                           _mx = - _mx, ## NOTE NEGATIVE 
                                           _my = _my,
                                           pulseEn = 2.51516e-2)))


            self.wfr.type = 'gsm'

        print("Coherent Modes Generated")
        self.log.info("Coherent Modes Generated")
        
        print("Testing Coherent Modes")
        self.log.info("Testing Coherent Modes")
        
        axis_x = np.linspace(-0.002, 0.002, 2048)
        axis_y = np.linspace(-0.002, 0.002, 2048)
        
        eField = np.zeros((2048,2048,1))
        
        for itr in range(len(self.wfr.modes)):
            eField += self.wfr.modes[itr].data.arrEhor[:,:,:,0] * self.wfr.eigenval[itr]
            [cohx, cohy] = coherence_log(eField, axis_x, axis_y)
            
            
            self.log.info("**************************************************")
            self.log.info("Addition of Mode {}".format(itr))
            self.log.info("Mode Weighting {}".format(self.wfr.eigenval[itr]))
            self.log.info("Mode Weighting {}".format(self.wfr.eigenfn[itr]))
            self.log.info("\n")
            self.log.info("Ix: {} m\nJx: {} m".format(cohx[2], cohx[3]))
            self.log.info("Iy: {} m\nJy: {} m".format(cohy[2], cohy[3]))
            self.log.info("x-beta: {} \n y-beta: {}".format(cohx[4], cohy[4]))

            res.append(cohx[2])
            res.append(cohy[2])
            res.append(cohx[3])
            res.append(cohy[3])
            res.append(cohx[4])
            res.append(cohy[4])

            results.append(res)
        
        return results

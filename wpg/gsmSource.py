import itertools

import numpy as np

from wpg.wavefront import Wavefront
from wpg.generators import build_gauss_wavefront_xy

from copy import copy, deepcopy
from operator import itemgetter

def addSlice(wfr, tmp_wfr):
    
    """
    adds the electric field of the temporary wavefront (and metadata) to the
    primary wfr structure
    
    :param wfr: primary wfr of gsmSource type
    :param wfr: temporary wfr of any wfr type
    """
    
    if wfr.params.type != 'Gaussian Schell Model Type Beam':
        print("wfr must be Gaussian Schell Model Type Beam")
    else:
        wfr.data.arrEhor = np.concatenate((wfr.data.arrEhor, tmp_wfr.data.arrEhor), axis = 2)
        wfr.data.arrEver = np.concatenate((wfr.data.arrEver, tmp_wfr.data.arrEver), axis = 2)
    
    wfr.params.Mesh.nSlices += 1
        
def gsmSource(wfr, mx, my, sigX = None, sigY = None, d2waist = None):
   
    """
    takes an existing wavefield structure and uses it to generate a gsm source
    w/ a coherence function defined by starikov and wolf
    """
    
    if sigX == None:
        sigX = (wfr.params.Mesh.xMax-wfr.params.Mesh.xMin)/10
    
    if sigY == None: 
        sigY = (wfr.params.Mesh.yMax-wfr.params.Mesh.yMin)/10
        
    wfr.params.type = 'Gaussian Schell Model Type Beam'
    
    #### GENERATE TRANSVERSE ELECTRIC MODE PAIRS
    wfr.params.eigenmodes = ordered_pairs(mx, my)
    
    #xMin, xMax, yMin, yMax = wfr.get_limits()
    for mx, my in wfr.params.eigenmodes:
        tmp_gsn = Wavefront(build_gauss_wavefront_xy(wfr.params.Mesh.nx, wfr.params.Mesh.ny,
                                                     wfr.params.photonEnergy/1000, *wfr.get_limits(), 
                                                     sigX, sigY, 
                                                     wfr.params.Mesh.zCoord,
                                                     _mx = mx, _my = my))
        addSlice(wfr,tmp_gsn)
 
    
    
    return wfr

class gsmSource():
    
    def __init__(self):
        """
        
        """
        self.wfr = Wavefront()
        self.eigenModes = []
        self.weightings = []

    
    
    def generatePairs(self, mx, my):
        """
        
        generates ordered pairs for use as transverse electric modes up (mx,my) 
        
        :param mx: highest horizontal transverse mode
        :param my: vertical transverse mode
        
        """
        pairs = []
        for x in range(mx+1):
            for y in range(my+1):
                pairs.append((x,y))
        
        self.eigenModes = pairs

         
    
    def generateEigenmodes(self, nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, d2Waist):
        """
        Generates eigenmodes of the partially coherent GSM beam
        
        :param nx: wfr horizontal dimensions in pixels (int)
        :param ny: wfr vertical dimensions in pixels (int)
        :param ekev: wfr energy in keV (float)
        :param xMin: wfr horizontal spatial minimum (float)
        :param xMax: wfr horizontal spatial maximum (float)
        :param yMin: wfr vertical spatial minimum (float)
        :param yMax: wfr vertical spatial maximum (float)
        :param sigX: fundamental gaussian horizontal FWHM (1st dev.) (float)
        :param sigY: fundamental gaussian vertical FWHM (1st dev.) (float)
        :param d2Waist: distance to waist - defines beam curvature (float)
        """
        
        
        
        for mx, my in self.eigenModes:
            
            if mx == 0 and my == 0:
                
                self.wfr = Wavefront(build_gauss_wavefront_xy(nx, ny,
                                                             ekev,
                                                             xMin, xMax, 
                                                             yMin, yMax, 
                                                             sigX, sigY, 
                                                             d2Waist,
                                                             _mx = mx, _my = my))
                
                self.wfr.params.type = 'Gaussian Schell Model Type Beam'
            else:
                
                tmp_gsn = Wavefront(build_gauss_wavefront_xy(nx, ny,
                                                             ekev,
                                                             xMin, xMax, 
                                                             yMin, yMax, 
                                                             sigX, sigY, 
                                                             d2Waist,
                                                             _mx = mx, _my = my))
                
                addSlice(self.wfr,tmp_gsn)
            
                del tmp_gsn
    
    def generateWeightings(self, sigma_I, sigma_mu, axisName = 'x'):
        """
        @author: gvr 
        @modded: twguest 25/03/20
        
        returns set of eigenvalues normalised to 1
        
        definitions follow Starikov and Wolf, 1982:
            
        beta =  p_mu/p_I is a measure of the "degree of global coherence" of the source.
       
        When beta >> 1, the source is effectively spatially coherent in the global sense and is then
        found to be well represented by a single mode. 
        
        When beta << 1, the source is effectively spatially incoherent in the global
        sense, and the number of modes needed to describe its behavior is larget
        
        :param sigma_I: 1D FWHM of intensity profile (float)
        :param sigma_mu: 1D Coherence Width (defined FWHM) (float)
        """
    
        a = 1 / (2 * sigma_I ** 2)
        b = 1 / (2 * sigma_mu ** 2)
        c = (a ** 2 + 2 * a * b) ** (1 / 2)
        
        weightfundamental = 1.0 ## Fundamental Mode Weight
        
        
        if axisName == 'x':
            M = max(self.eigenModes,key=itemgetter(0))[0] # gets max x mode
        elif axisName == 'y':
            M = max(self.eigenModes,key=itemgetter(1))[1] # gets max y mode
        
        
        self.weightings = [weightfundamental * (b / (a + b + c)) ** n for n in range(M+1)]
        self.weightings /= np.sum(self.weightings)
        
    
    def generateWeights2D(self, sigma_Ix, sigma_mux, sigma_Iy, sigma_muy):
        """
        @author: twguest 
        @modded: twguest 25/03/20
        
        returns set of eigenvalue normalised to 1
        
        under the assumption that the transverse components of the MCF are seperable
        definitions follow Starikov and Wolf, 1982:
            
        beta =  p_mu/p_I is a measure of the "degree of global coherence" of the source.
       
        When beta >> 1, the source is effectively spatially coherent in the global sense and is then
        found to be well represented by a single mode. 
        
        When beta << 1, the source is effectively spatially incoherent in the global
        sense, and the number of modes needed to describe its behavior is larget
        
        :param sigma_Ix: 1D FWHM of horizontal intensity profile (float)
        :param sigma_mux: 1D horizontal Coherence Width (defined FWHM) (float)
        
        :param sigma_Iy: 1D FWHM of vertical intensity profile (float)
        :param sigma_muy: 1D vertical Coherence Width (defined FWHM) (float)
        """
        
        a = 1 / (2 * sigma_I ** 2)
        b = 1 / (2 * sigma_mu ** 2)
        c = (a ** 2 + 2 * a * b) ** (1 / 2)
        
        weightfundamental = 1.0 ## Fundamental Mode Weight
        
        

        Mx = max(a,key=itemgetter(0))[0] # gets max x mode
        My = max(a,key=itemgetter(1))[1] # gets max y mode
    
       
        self.weightings = [weightfundamental * (b / (a + b + c)) ** n for n in range(M)]
        self.weightings /= np.sum(self.weightings)

    
    def replaceWavefront(self, new_wfr):
        """
        replaces current wfr with new_wfr, useful after propagation
        """
        self.wfr = new_wfr
        
    def collapseModes(self):
        """
        (this should be checked)
        
        :returns cmplx: scaled cmplx incoherent sum of modes valued array
        """
        wfield = []
        
        cmplx = self.wfr.toComplex()
        print(cmplx.shape)
        for pol in cmplx:
            print(pol.shape)
            for itr in range(pol.shape[2]):
                pol[:,:,itr] *= self.weightings[itr]
            
            wfield.append(np.sum(pol, axis = 2))
            
        return wfield      
                    

    def get_wavefront(self):
        """
        :returns self.wfr: Wavefront object of GSM source
        """
        return self.wfr
    
    def get_fundamental(self):
        """  
        :returns fnd: [Wavefront Obj.] 
        fundamental mode of gaussian schell model source TEM_{00}.
        """
        fnd = deepcopy(self.wfr)
        
        fnd.params.type = 'fundamental mode'
        
        fnd.data.wfr.arrEhor = fnd.data.wfr.arrEhor[:,:,0,:] 
        fnd.data.wfr.arrEver = fnd.data.wfr.arrEver[:,:,0,:]
        
        return fnd
    
def mcf(E):
    """
    calculate mutual coherence function of the wavefield relative to the central
    point
    
    :param E: complex wavefield [x,y,slices]
    
    :returns mcf: mutual coherence function    
    """
    
    max_index = np.unravel_index(E.argmax(), E.shape)
    
    Emax = np.max(E)
    
    MCF = np.mean(Emax*np.conj(E), axis = 2)
    
    return MCF

def cdoc(E):
    """
    calculate complex degree of coherence
    
    :param E: complex wavefield [x,y,slices]
    
    returns: cdoc - complex degree of coherence
    """
    
    m = mcf(E)
    Emax = np.max(E)
    
    cdoc = m/(np.sqrt(Emax**2)*np.sqrt(np.mean(E**2,axis = 2)))
    
    return cdoc

def csd(E):
    
    

    MCF = mcf(E)
    csd = (np.fft.fft2(MCF))
    
    
    
    S1 =  np.fft.fft2(np.mean(abs(E)**2, axis = 2))
    S2 =  (np.mean(abs(np.max(E))**2))
    
    sdoc = csd/(np.sqrt(S2)*np.sqrt(S1))
    return sdoc


def tdoc(E):
    nx, ny, slc = E.shape
    k = nx*ny
    D = np.array(E).reshape((slc, k), order='F').T
    DTD = np.dot(D.T.conjugate(), D)
    res = np.diag(np.dot(DTD, DTD)).sum() / np.diag(DTD).sum()**2
    return res.real

if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    
    gsmSource = gsmSource()
    gsmSource.generatePairs(3, 0)
    gsmSource.generateEigenmodes(256, 256, 9.2, -5e-07, 5e-07, -5e-07, 5e-07, 1e-07, 1e-07, 1e-09)
    gsmSource.generateWeightings(1e-07, 1e-08)
    E = gsmSource.collapseModes()
    
    plt.imshow(np.abs(E[0]**2))
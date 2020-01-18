import numpy as np

from scipy.signal import correlate2d as correlate
from scipy.signal import argrelextrema as extrema
from sklearn.preprocessing import minmax_scale as norm
from extensions.twpg_uti_wf import fixPhase
from copy import copy

def tdoc(E):
    
    repeats, binsx, binsz = U.shape
    k = binsx * binsz
    D = np.array(U).reshape((repeats, k), order='F').T
    DTD = np.dot(D.T.conjugate(), D)
    res = np.diag(np.dot(DTD, DTD)).sum() / np.diag(DTD).sum()**2
    return res.real

def beta(mu_, sig):
    beta = mu_/sig
    return beta

def intensity_width(II, axis, thr = 0.88):
    
    print("Calculating Intensity Profile Width")
    var = np.sum((I - np.mean(II))**2 for I in II)/np.sum(II)
    std = np.sqrt(var)
    
    sig = extrema(II, np.less)[0]
    sig = sig[(axis[sig] > len(II)//2) & (II[sig][:,0] < np.mean(II) + std*2)]
    
    sig = axis[sig[0]] - axis[len(II)//2]
    
    return sig

def coherence_len(E, axis, thr = 0.88):
    
    j12 = coherence_1d(E)
    
    print("Calculating Coherence Length")
    mu_ = j12[len(j12)//2:]
    mu_ = extrema(j12, np.less)[0]
    
    mu_ = mu_[(mu_ > len(j12)//2) & (j12[mu_][0]< thr)]
    print(mu_)
    mu_ = axis[mu_[0]] - axis[len(j12)//2]
    
    return j12, mu_

def coherence_angle(E):
    phase = copy(E.imag)
    mid = copy(E.imag)
    mid[:] = E.imag[len(E)//2]
    diff = []
    for itr in range(len(E)):
        d_ = np.cos(abs(mid[itr]-phase[itr]))
        diff.append(d_)
    print(np.abs(diff[len(E)//2:len(E)//2+50]))


def coherence_xrt(E):
    print("Trying XRT Method")
    j12 = []
    J = np.dot(E ,E.T.conjugate())

    II = np.diag(J)
    J /= II**0.5 * II[:,np.newaxis]**0.5
    Jd = np.diag(np.fliplr(J))
    print(Jd)
    return Jd
def coherence_1d(E):
    
    print("Calculating 1D DoC")
    j12 = []
    
    for itr in range(len(E)):
        J12 = np.dot(E[itr], E[len(E)//2].conjugate())
        J11 = np.dot(E[itr], E[itr].conjugate())
        J22 = np.dot(E[len(E)//2], E[len(E)//2].conjugate())
        j12.append(J12/(np.sqrt(J11 * J22)))
        
            
    j12 = np.abs(j12)*np.cos(np.imag(j12))
    
    return j12

if __name__ == '__main__':
    
    from extensions.twpg_wavefront import Wavefront
    
    
    for itr in range(1):
        
        if itr == 0:
            wfr = Wavefront()
            wfr.load_hdf5(r"/nfs/data/users/twg/gsmProp/atSource/wfr_mode_{}.hdf5".format(itr))
            eField = copy(wfr.data.arrEhor)
            eField[:,:,0,0] = norm(np.nan_to_num(eField[:,:,0,0]))
        else:
            wfr = Wavefront()
            wfr.load_hdf5(r"/nfs/data/users/twg/gsmProp/atSource/wfr_mode_{}.hdf5".format(itr))
            eField[:,:,0,0] += norm(np.nan_to_num(wfr.data.arrEhor[:,:,0,0]))
            eField[:,:,0,1] += copy(wfr.data.arrEhor[:,:,0,1]) 
    
    
    
    x,y,z,c = eField.shape
    fixPhase(eField)

    E = np.zeros((x,y,1)).astype(complex)
    E[:,:, 0] += eField[:,:,0,0]
    E[:,:, 0] += eField[:,:,0,1]*1j
            
    ### MAKE CUTS
    E = E[:,y//2,:] #xcut
    
    print("---------------------------")
    coherence_angle(E)

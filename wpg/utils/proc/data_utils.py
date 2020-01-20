# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import hilbert, chirp
from extensions.wpg_uti_wf import *
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from sklearn.preprocessing import minmax_scale
from scipy import optimize
from pylab import *
from scipy.signal import argrelextrema

def lineprof(wfr, axis = 'x', buffer = 5):
    """
    Returns a line profile across the centre of the chosen axis of a wavefront
    
    :param wfr: wfr
    :param axis: axis of choice
    """
    
    if axis == 'x':    
        mid = wfr.get_intensity().shape[0]//2
        profile = wfr.get_intensity()[mid,:,0]
        
        for i in range(1,buffer+1):
            profile += wfr.get_intensity()[mid+buffer,:,0]
            profile += wfr.get_intensity()[mid-buffer,:,0]
            profile /= 2*buffer + 1

    elif axis == 'y':
        mid = wfr.get_intensity().shape[1]//2
        profile = wfr.get_intensity()[:,mid,0]
    
        for i in range(1,buffer+1):
            profile += wfr.get_intensity()[:,mid+buffer,0]
            profile += wfr.get_intensity()[:,mid-buffer,0]
            profile /= 2*buffer + 1
    else:
        print("Incorrect Axis: Should be 'x' or 'y'")
    
    return profile


def model_lineprof(wfr, axis = 'x'):
    
    if axis == 'x':
        profile = lineprof(wfr, 'x')
        profile = minmax_scale(profile)
        x = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, wfr.params.Mesh.nx)
        x2 = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, 1500)

    elif axis == 'y':
        profile = lineprof(wfr, 'y')
        profile = minmax_scale(profile)
        x = np.linspace(wfr.params.Mesh.yMin, wfr.params.Mesh.yMax, wfr.params.Mesh.yx)
    else:
        print("Incorrect Axis: Should be 'x' or 'y'")

    env = np.abs(hilbert(profile, N = wfr.params.Mesh.nx))

    envmm = argrelextrema(env, np.less)[0]
    print(len(envmm))
# =============================================================================
# 
# 
#     mini = min(envmm[np.nonzero(envmm)])
#     plus = env[len(env)//2 + mini]
#     minus = env[len(env)//2 - mini]
#     
#     if plus <= minus:
#         minimum = plus
#     elif minus < plus:
#         minimum = minus
#         
#     ### get max
#     maximum = env[0]
#         
# =============================================================================
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, len(env)), env)
    fig.savefig("test")
    return profile, env

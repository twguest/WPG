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
import seaborn as sns    
from extensions.twpg_wavefront import Wavefront
from extensions.twpg_uti_wf import fixPhase

font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 14}

plt.rc('font', **font)

from math import trunc
#plt.style.use('seaborn-white')
import cplot
import numexpr as ne
from numpy import pi


def plotWavefront(wfr):
    sns.set()
        
    wfr.type = 'se'
    wfr.electric_field()
    wfr.intensity()
    fixPhase(wfr.eField)
    
    ax = sns.heatmap(wfr.II[:,:,0])
    
    return ax
    

if __name__ == "__main__":
    print("Testing Plotting Utility")
    
    wfr = Wavefront()
    wfr.load_hdf5(r"/nfs/data/users/twg/gsmProp/9_SRWLOptD/wfr_mode_0.hdf5")
    
    
    ax = plotWavefront(wfr)
    fig = ax.get_figure()
    fig.savefig(r"./extensions/out/preES")

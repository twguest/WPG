#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:21:46 2019

@author: alaleh
"""


import h5py
from wavefront import Wavefront
import numexpr as ne
import csv
from extensions.gvrutils import propagationParameters, Struct
from extensions.beamlineExtAl2 import bl

from wpg.srwlib import SRWLOptC, SRWLOptD


path = '/home/alaleh/Desktop/out/2019-07-16/64dd8f7e-a759-11e9-a3a0-1f9506ef013f/wfPostOSA.h5'  #/home/alaleh/Desktop/out/2019-07-10/30Slices_63D_312nmt/wfPostOSA.h5'
px = 2.205e-8 # nm
ny, nx = 128*8,128*8  # output dimensions
outPath = path[:-2] + '1024.h5'

wfr=Wavefront()
wfr.load_hdf5(path)


Nx, Ny = wfr.params.Mesh.nx, wfr.params.Mesh.ny


Resc = Struct()
Resc.driftLength=0.0
Resc.pp =propagationParameters(SemiAnalyt = 0.0, RangeX=nx/Nx, RangeY=ny/Nx, ResolutionX = 1., ResolutionY = 1. )  #Res=147 to get 10nm, Res=73.5 to get 20nm #183.74 to get 8nm # 168 to get 8.75nm
Resc.description='Rescale only'


elResc = bl(SRWLOptC(_arOpt=[SRWLOptD(Resc.driftLength)], 
                    _arProp=[Resc.pp]), 
                    description=Resc.description)

elResc.propagate(wfr=wfr,
                  plot=True,
                  outFilePath=None, 
                  propagationOverride=None)


re = wfr.get_real_part()
im = wfr.get_imag_part()
cplx = ne.evaluate("complex(re,im)")

f = h5py.File(outPath,'w')
f['probe'] = cplx
f.close

f = h5py.File(outPath,'r')
cplx = f['probe'] 
f.close





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 21:09:41 2019

@author: gvanriessen
"""


#    PTYCHO=False
#    if PTYCHO==True:
#
#        # object specification.  3 um thick PMMA @ 460 eV
#        obj=Struct()
#        #######obj.file_path,
#        obj.resolution=10.e-9,
#        obj.thickness=3.e-6,
#        obj.delta=0.00112663  # C5H8O2 Density=1.19, Energy=460 eV, Beta = 0.000272084551
#        obj.atten_len= 0.788312
#        obj.xc = 0.0
#        obj.yc = 0.0
#        obj.area=None
#        obj.rotate_angle=None
#        obj.rotate_reshape=False
#        obj.cutoff_background_noise=0
#        obj.background_color=0
#        obj.tile=None
#        obj.invert=False
#        obj.is_save_images = False
#
#        # propagate to object plane
#        drift = blDrift(600.e-6, RangeX=1.0, RangeY=1.0, ResolutionX=2.0, ResolutionY=2.0)
#        drift.propagate(wfr)
#
#        #position list
#        scanRange=10.e-6
#        positions = snakeScan(centre=(0,0), range=(scanRange,scanRange))
#        positions = [[0,0]] #(test with one position:
#
#
#        for pos in positions:
#            wf=wfr
#
#            shiftX, shiftY = pos
#
#            #define object transmission function & it's position
#            obj.shift_x, obj.shift_y = pos
#            tr = objOpt (obj)
#
#            el.append(tr)

#    el, pp = [], []
#
#    #drift 5m
#    el.append(SRWLOptD(5))
#    pp.append(propagationParameters(RangeX=2.0, RangeY=2.0, ResolutionX=2.0, ResolutionY=2.0))
#
#
#    # Aperture
#    el.append(SRWLOptA("r", "a", 100e-6, 100e-6, 0.0, 0.0))
#    pp.append(propagationParameters(RangeX=0.1, RangeY=0.1, ResolutionX=10.0, ResolutionY=10.0))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:41:58 2019

@author: Grant van Riessen

"""

#!/usr/bin/env python
import os, subprocess, re, uuid, sys
import multiprocessing as mp
import numpy as np

def run(job): 
    cmd, paramFile, S, cwd = job  # unpack arguments    
    
    print ('Running {} for parameters in {}, scaling mask RMS by {}.  Working directory = {}'.format(cmd, paramFile, S, cwd))        
    
    subprocess.check_call([ cmd + ' ' + paramFile +  ' ' + S ],cwd=cwd, shell=True)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: run(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def main():
    
    # working directory
    workdir =r'/user/home/WPG/'
    
    # path to file containing parameters, relative to workdir
    parameterPath = 'extensions.parameters.BLParams_Kino'
    
    # cmd to execute
    script = 'source activate py36; cd /user/home/WPG; ipython /user/home/WPG/extensions/blKinoform2.py '
    
    jobs = []
    for i in np.linspace(0,2,10): # i is roughness scaling factor
        jobs.append([script, parameterPath, str(i), workdir])

    # start processes
    pool = mp.Pool() # use all available CPUs
    pool.map(safe_run, jobs)

if __name__=="__main__":
    mp.freeze_support() # optional if the program is not frozen
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:41:58 2019

@author: Grant van Riessen

"""

#!/usr/bin/env python
import os, subprocess, re
import multiprocessing as mp

def run(cmd_param_cwd): 
    cmd, paramFile, cwd = cmd_param_cwd # unpack arguments    
    print ('Running {} for parameters in {}.  Working directory = {}'.format(cmd, paramFile, cwd))        
    
    os.chdir(cwd)
    subprocess.check_call([cmd, paramFile],cwd=cwd)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: run(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def main():
    
    
    # working directory
    workdir =r'/opt/wpg/wpg/'
    
    # directory containing parameters, relative to workdir
    parameterDir = 'extensions/parameters/batch0/'
    
    # script to execute
    script = r'/opt/wpg/wpg/extensions/blLTZP0.py'
    
    files = []
    for f in os.listdir(os.path.join(workdir, parameterDir)):
        if f.endswith('.py'):
           arg = re.sub("/", ".", parameterDir) + f[:-3]
           files.append([script, arg, workdir])

    # start processes
    pool = mp.Pool() # use all available CPUs
    pool.map(safe_run, files)

if __name__=="__main__":
    mp.freeze_support() # optional if the program is not frozen
    main()
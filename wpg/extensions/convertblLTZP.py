#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:22:59 2019

@author: gvanriessen
"""

import IPython.nbformat.current as nbf
nb = nbf.read(open('blLTZPeff.py', 'r'), 'py')
nbf.write(nb, open('blLTZPeff.ipynb', 'w'), 'ipynb')
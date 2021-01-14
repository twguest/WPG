# -*- coding: utf-8 -*-
"""
This module contains utils and classes for X-ray wavefront propagation.

The most propagation methods based on SRW library.

.. module:: wpg
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# This deprecate warnings form SRWLib visualization module

# import warnings
# warnings.filterwarnings("ignore")
# import srwlpy as srwlpy
# import srwlib as srwlib
# warnings.resetwarnings()

# fix segmentation fault using fftw from numpy mkl
# from . import srwlpy

# Create aliases for simple importing

from .wavefront import Wavefront
from .beamline import Beamline
from .glossary import print_glossary


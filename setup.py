#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:38:30 2020

@author: twguest
"""

from setuptools import Extension
from distutils.core import setup, Extension

from distutils.command.install import install as DistutilsInstall
from os import system
from distutils.core import setup

class MyInstall(DistutilsInstall):
    def run(self):
        
        system("make all")
        DistutilsInstall.run(self)

        

setup(
    name='WPG',
    version='0.1.0',
    packages=['wpg'],
    license='MIT',
    long_description=open('README.md').read()
    cmdclass={'install': MyInstall},
)

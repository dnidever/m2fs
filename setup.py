#!/usr/bin/env python

from distutils.core import setup

setup(name='m2fs',
      version='1.0',
      description='M2FS Data Reduction Pipeline',
      author='Matt Walker, David Nidever',
      author_email='mgwalker21@gmail.com',
      url='https://github.com/dnidever/m2fs',
      packages=['m2fs'],
      scripts=['bin/m2fs_zero','bin/m2fs_dark','bin/m2fs_reduce','bin/m2fs_apertures'],
      requires=['numpy','astropy','scipy','specutils','ccdproc']
)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#################################################################################
#                                                                               # 
# TDOC - Thermodynamic Data Off-line Calculator                                 #
#                                                                               #
# Copyright (c) 2022 Center for Combustion Dynamics, Sichuan University.        #
#                                                                               #
# Permission is hereby granted, free of charge, to any person obtaining a copy  # 
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      # 
#                                                                               #
# The above copyright notice and this permission notice shall be included in    #
# all copies or substantial portions of the Software.                           #
#                                                                               #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE #
# SOFTWARE.                                                                     #
#                                                                               #
#################################################################################

from os import system
from setuptools import setup



try:
      import rdkit, ase, scipy, numpy, openbabel
except:
    system('conda install --yes -c conda-forge --file requirements.txt')
finally:
      import numpy


setup(name = 'tdoc',
      version = '1.0',
      description = 'Thermodynamic data offline calculator',
      packages = ['tdoc'],
      scripts = ['tdoc.py'],
      license = 'MIT',
      platforms = ['Windows', 'Linux'],
      url = 'https://github.com/lixiangyuangroup/TDOC',
      author_email = 'renhs@scu.edu.cn',
      author = 'Huajie Xu, Zihan Xu, Lu Liu, Zerong Li, Quan Zhu, Haisheng Ren',
      include_dirs=['.', numpy.get_include()],
      )

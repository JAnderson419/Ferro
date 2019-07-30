#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from os.path import join, dirname, realpath
from context import models as lf
from context import data as hd

plt.close('all')
discretecapdir = join(dirname(realpath(__file__)), 'testData', '2pt2nF_discrete_freqs')
dcfiles = hd.dir_read(discretecapdir)
data = hd.list_read(dcfiles)
testfilm = lf.LandauFilm()   
cde = testfilm.c_calc(data, plot=1)
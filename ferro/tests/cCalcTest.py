#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')
discretecapdir = r'.\testData\2pt2nF_discrete_freqs'
dcfiles = hd.dirRead(discretecapdir)
data = hd.listRead(dcfiles)
testfilm = lf.LandauFilm()   
cde = testfilm.c_calc(data, plot=1)
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
testfilm = lf.LandauFilm()   
cde = testfilm.cCalc(dcfiles, plot=1)
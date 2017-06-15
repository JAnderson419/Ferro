#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')
discretecapdir = r'.\testData\MDCref\mdc100pf'
files = hd.dirRead(discretecapdir)
data = hd.listRead(files)
testfilm = lf.LandauFilm()
cde = testfilm.cCalc(data, plot=1)
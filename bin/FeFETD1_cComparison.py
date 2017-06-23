#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')
freqdirs = [r"..\ferro\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_freqs",
            r'..\ferro\tests\testData\FeFETD1\MFS+\FeFETD1_die68_MFS+_100_10x10_freqs'
            ]
for f in freqdirs:
    dcfiles = hd.dirRead(f)
    data = hd.listRead(dcfiles)
    testfilm = lf.LandauFilm()   
    cde = testfilm.cCalc(data, plot=1)
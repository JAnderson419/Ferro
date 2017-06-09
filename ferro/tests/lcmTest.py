#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson

Tests fftPlot and bandstopFilter functions of HysteresisData
"""

import matplotlib.pyplot as plt
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')

tempdir = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_temps"
tempfiles = hd.dirRead(tempdir)
templkgdir = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_tempslkg"
templkgfiles = hd.dirRead(templkgdir)

datafile = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_temps\H9 die (9,4) S3 79C 100Hz 3V 1Average Table3.tsv"
lkgfile = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_tempslkg\H9 die (9,4) S3 79C 2s step Table2.tsv"

data = hd.HysteresisData()
data.tsvRead(datafile)

ldata = hd.LeakageData()
ldata.lcmRead(lkgfile)

data.hystPlot()
ldata.lcmFit()
ldata.lcmPlot()
compData = data.leakageCompensation(ldata)
compData.hystPlot()
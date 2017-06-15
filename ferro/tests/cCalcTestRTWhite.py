#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')


RTfreqDir = r".\testData\RT WhiteA\RTWhiteAFreq"
RTfreqFiles = hd.dirRead(RTfreqDir)
RTfreqData = hd.listRead(RTfreqFiles)
RTfreq100hz = r".\testData\RT WhiteA\RTWhiteAFreq\RT WhiteA 100Hz 8V 1Average Table1.tsv"

RT100data = hd.HysteresisData()
RT100data.tsvRead(RTfreq100hz)
RT100data.hystPlot

RTWhiteFilm = lf.LandauSimple(thickness = 255E-7, area=1E-4)
RTWhiteFilm.c = RTWhiteFilm.cCalc(RTfreqData, plot = 1)
RT100compensated, RTWhiteFilm.pr = RTWhiteFilm.cCompensation(RT100data, plot = 1)


# Following code plots a series of diff freq hystdata files on same plot

hystData = []
legend = []
for f in RTfreqFiles:
    data = hd.HysteresisData()
    data.tsvRead(f)
    data.dvdtPlot() # plots dvdt for analysis - unrelated to freq hystPlot
    hystData.append(data)
    legend.append(int(data.freq))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.freq))  

legend = [str(x)+' Hz' for x in legend]      
hd.hystPlot(hystData, legend)
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
RTfreqFiles = hd.dir_read(RTfreqDir)
RTfreqData = hd.list_read(RTfreqFiles)
RTfreq100hz = r".\testData\RT WhiteA\RTWhiteAFreq\RT WhiteA 100Hz 8V 1Average Table1.tsv"

RT100data = hd.HysteresisData()
RT100data.tsv_read(RTfreq100hz)
RT100data.hyst_plot

RTWhiteFilm = lf.LandauSimple(thickness = 255E-7, area=1E-4)
RTWhiteFilm.c = RTWhiteFilm.c_calc(RTfreqData, plot = 1)
RT100compensated, RTWhiteFilm.pr = RTWhiteFilm.c_compensation(RT100data, plot = 1)


# Following code plots a series of diff freq hystdata files on same plot

hystData = []
legend = []
for f in RTfreqFiles:
    data = hd.HysteresisData()
    data.tsv_read(f)
    data.dvdt_plot() # plots dvdt for analysis - unrelated to freq hyst_plot
    hystData.append(data)
    legend.append(int(data.freq))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.freq))  

legend = [str(x)+' Hz' for x in legend]      
hd.hyst_plot(hystData, legend)
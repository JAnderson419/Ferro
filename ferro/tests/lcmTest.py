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



RTfreq1000hz = r".\testData\RT WhiteA\RT WhiteA 1Hz 8V 1Average Table3.tsv"

RTdata = hd.HysteresisData()
RTdata.tsvRead(RTfreq1000hz)
RTdata.hystPlot()

RTdata.fftPlot(RTdata.current)
RTdata.current = RTdata.bandstopFilter(RTdata.current, plot = True)
RTdata.polarization = RTdata.bandstopFilter(RTdata.polarization)
RTdata.fftPlot(RTdata.current)
#RTdata.lpFilter(RTdata.polarization)
RTdata.hystPlot()


# Following code plots a series of diff temp hystdata files on same plot

hystData = []
legend = []
for f in tempfiles:
    data = hd.HysteresisData()
    data.tsvRead(f)
    hystData.append(data)
    legend.append(int(data.temp))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.temp))

legend = [str(x)+' K' for x in legend]  
hd.hystPlot(hystData, legend)
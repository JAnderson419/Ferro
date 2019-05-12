#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson

Tests fft_plot and bandstop_filter functions of HysteresisData
"""

import matplotlib.pyplot as plt
from os.path import join, dirname, realpath
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')



RTfreq1000hz = join(dirname(realpath(__file__)), 'testData', 'RT WhiteA', 'RT WhiteA 1Hz 8V 1Average Table3.tsv')

RTdata = hd.HysteresisData()
RTdata.tsv_read(RTfreq1000hz)
RTdata.hyst_plot()

RTdata.fft_plot(RTdata.current)
RTdata.current = RTdata.bandstop_filter(RTdata.current, plot=True)
RTdata.polarization = RTdata.bandstop_filter(RTdata.polarization)
RTdata.fft_plot(RTdata.current)
#RTdata.lpFilter(RTdata.polarization)
RTdata.hyst_plot()
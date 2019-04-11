#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')

freqdir = r".\testData\hfo2_MFM\H9_x9y4_1e4_freq"
freqfiles = hd.dirRead(freqdir)
freqdata = hd.listRead(freqfiles)
hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)

hfo2.c = hfo2.c_calc(freqdata, plot=1)
hfo2.rho_calc(freqdata)

tempdir = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_temps"
tempfiles = hd.dirRead(tempdir)
tempdata = hd.listRead(tempfiles)
templkgdir = r".\testData\hfo2_MFM\H9_x9y4_1e4_S3_tempslkg"
templkgfiles = hd.dirRead(templkgdir)

hfo2.a0 = hfo2.a0_calc(tempdata)


# Following code plots a series of diff freq hystdata files on same plot

hystData = []
legend = []
for f in freqfiles:
    data = hd.HysteresisData()
    data.tsv_read(f)
#    data.dvdt_plot() # plots dvdt for analysis - unrelated to freq hyst_plot
    hystData.append(data)
    legend.append(int(data.freq))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.freq))

legend = [str(x)+' Hz' for x in legend]  
hd.hyst_plot(hystData, legend)

# Following code plots a series of diff temp hystdata files on same plot

hystData = []
legend = []
for f in tempfiles:
    data = hd.HysteresisData()
    data.tsv_read(f)
    hystData.append(data)
    legend.append(int(data.temp))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.temp))

legend = [str(x)+' K' for x in legend]  
hd.hyst_plot(hystData, legend)

# Following code plots a series of diff temp leakagedata files on same plot

leakageData = []
legend = []
for f in templkgfiles:
    data = hd.LeakageData()
    data.lcm_read(f)
    leakageData.append(data)
    legend.append(int(data.temp))

legend = sorted(legend)
leakageData = sorted(leakageData, key=lambda data: int(data.temp))
legend = [str(x)+' K' for x in legend]  
hd.lcm_plot(leakageData, legend)
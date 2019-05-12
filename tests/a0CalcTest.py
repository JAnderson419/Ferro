#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from os.path import join
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')


sampledir = join('.', 'testData', 'hfo2_MFM')
freqdir = join(sampledir, 'H9_x9y4_1e4_freq')
freqfiles = hd.dir_read(freqdir)
freqdata = hd.list_read(freqfiles)
hfo2 = lf.LandauFull(thickness=13E-7, area=6579E-8)

hfo2.c = hfo2.c_calc(freqdata, plot=1)
hfo2.rho_calc(freqdata)

tempdir = join(sampledir, 'H9_x9y4_1e4_S3_temps')
tempfiles = hd.dir_read(tempdir)
tempdata = hd.list_read(tempfiles)
templkgdir = join(sampledir, 'H9_x9y4_1e4_S3_tempslkg')
templkgfiles = hd.dir_read(templkgdir)

hfo2.a0 = hfo2.a0_calc(tempdata)


# Following code plots a series of diff freq hystdata files on same plot

hyst_data = []
legend = []
for f in freqfiles:
    data = hd.HysteresisData()
    data.tsv_read(f)
#    data.dvdt_plot() # plots dvdt for analysis - unrelated to freq hyst_plot
    hyst_data.append(data)
    legend.append(int(data.freq))

legend = sorted(legend)
hyst_data = sorted(hyst_data, key=lambda data: int(data.freq))

legend = [str(x)+' Hz' for x in legend]  
hd.hyst_plot(hyst_data, legend)

# Following code plots a series of diff temp hystdata files on same plot

hyst_data = []
legend = []
for f in tempfiles:
    data = hd.HysteresisData()
    data.tsv_read(f)
    hyst_data.append(data)
    legend.append(int(data.temp))

legend = sorted(legend)
hyst_data = sorted(hyst_data, key=lambda data: int(data.temp))

legend = [str(x)+' K' for x in legend]  
hd.hyst_plot(hyst_data, legend)

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
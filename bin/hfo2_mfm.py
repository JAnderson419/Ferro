#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from context import models as lf
from context import data as hd
from os.path import join, dirname, realpath


plt.close('all')
testdatadir = join(dirname(dirname(realpath(__file__))), "tests", "testData")

freqdir = join(testdatadir, r"hfo2_MFM", "H9_x9y4_1e4_freq")

tempdir = join(testdatadir, r"hfo2_MFM", "H9_x9y4_1e4_S3_temps")
templkgdir = join(testdatadir, r"hfo2_MFM", "H9_x9y4_1e4_S3_tempslkg")
forcFile = join(testdatadir, r"hfo2_MFM", "H9_x9y4_1e4_forc",
                "H9 die (9,4) 0Hz 4V 1Average Table1.tsv")

templkgfiles = hd.dir_read(templkgdir)

tempfiles = hd.dir_read(tempdir)
tempData = hd.list_read(tempfiles, templkgfiles)

freqfiles = hd.dir_read(freqdir)
freqData = hd.list_read(freqfiles)
hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)
cCompData = freqData[0]

hfo2.c = hfo2.c_calc(freqData, plot=1)
compensatedData, hfo2.pr = hfo2.c_compensation(cCompData)
compensatedData.hyst_plot(plot_e=True)
hfo2.rho_calc(freqData)

hfo2.a0 = hfo2.a0_calc(tempData)

freqDataLkgComp = hd.list_read(freqfiles, templkgfiles)
cCompDataLkgComp = freqDataLkgComp[0]
hd.hyst_plot([cCompData, cCompDataLkgComp],
             ["With Leakage", "Without Leakage"], plot_e=1)

### FORC Calculation


hfo2_forc = hd.HysteresisData(area=6579E-8, thickness=13E-7)
hfo2_forc.tsv_read(forcFile)
hfo2_forc.hyst_plot(plot_e=1)
e, er, probs = hfo2_forc.forc_calc(plot = False)
    
domains = hfo2.domain_gen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-4.5E6,4.5E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
hfo2.calc_efe_preisach(esweep, domains, plot=1)

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
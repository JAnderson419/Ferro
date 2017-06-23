#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
from context import LandauFilm as lf
from context import HysteresisData as hd



plt.close('all')

freqdir = r"..\ferro\tests\testData\hfo2_MFM\H9_x9y4_1e4_freq"
tempdir = r"..\ferro\tests\testData\hfo2_MFM\H9_x9y4_1e4_S3_temps"
templkgdir = r"..\ferro\tests\testData\hfo2_MFM\H9_x9y4_1e4_S3_tempslkg"
forcFile = r"..\ferro\tests\testData\hfo2_MFM\H9_x9y4_1e4_forc\H9 die (9,4) 0Hz 4V 1Average Table1.tsv"

templkgfiles = hd.dirRead(templkgdir)

tempfiles = hd.dirRead(tempdir)
tempData = hd.listRead(tempfiles, templkgfiles)

freqfiles = hd.dirRead(freqdir)
freqData = hd.listRead(freqfiles)
hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)
cCompData = freqData[0]

hfo2.c = hfo2.cCalc(freqData, plot=1)
compensatedData, hfo2.pr = hfo2.cCompensation(cCompData)
compensatedData.hystPlot(plotE=True)
hfo2.rhoCalc(freqData)

hfo2.a0 = hfo2.a0Calc(tempData)

freqDataLkgComp = hd.listRead(freqfiles, templkgfiles)
cCompDataLkgComp = freqDataLkgComp[0]
hd.hystPlot([cCompData,cCompDataLkgComp],
            ["With Leakage","Without Leakage"],plotE=1)

### FORC Calculation


hfo2_forc = hd.HysteresisData(area=6579E-8, thickness=13E-7)
hfo2_forc.tsvRead(forcFile)
hfo2_forc.hystPlot(plotE=1)
e, er, probs = hfo2_forc.forcCalc(plot = False)
    
domains = hfo2.domainGen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-4.5E6,4.5E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
hfo2.calcEfePreisach(esweep, domains, plot=1)

# Following code plots a series of diff freq hystdata files on same plot

hystData = []
legend = []
for f in freqfiles:
    data = hd.HysteresisData()
    data.tsvRead(f)
#    data.dvdtPlot() # plots dvdt for analysis - unrelated to freq hystPlot
    hystData.append(data)
    legend.append(int(data.freq))

legend = sorted(legend)
hystData = sorted(hystData, key=lambda data: int(data.freq))

legend = [str(x)+' Hz' for x in legend]  
hd.hystPlot(hystData, legend)

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

# Following code plots a series of diff temp leakagedata files on same plot

leakageData = []
legend = []
for f in templkgfiles:
    data = hd.LeakageData()
    data.lcmRead(f)
    leakageData.append(data)
    legend.append(int(data.temp))

legend = sorted(legend)
leakageData = sorted(leakageData, key=lambda data: int(data.temp))
legend = [str(x)+' K' for x in legend]  
hd.lcmPlot(leakageData, legend)
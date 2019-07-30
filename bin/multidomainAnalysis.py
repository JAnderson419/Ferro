#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
import numpy as np
from context import models as lf
from context import data as hd

plt.close('all')

leakageComp = False
device = 0

### Radiant Technologies White B ###
if device == 0:
    freqdir = r"..\tests\testData\RTWhiteB\RTWhiteB_freqs"
    tempdir = None
    templkgdir = r"..\tests\testData\RTWhiteB\RTWHITEB_lkg"
    forcFile = r"..\tests\testData\RTWhiteB\RTWhiteB_FORC\RTWhiteB 0Hz 5V 1Average Table1.tsv"
    t = 255E-7 
    a = 1E-4 # mask defined area that was used in measurement 
    aReal = 1E-4 # includes effect of undercut during M1 etch

### Radiant Technologies White A ###
if device == 1:
    freqdir = r"..\tests\testData\RT WhiteA\RTWhiteAFreq"
    tempdir = None
    templkgdir = r"..\tests\testData\RTWhiteB\RTWHITEB_lkg"
    forcFile = r"..\tests\testData\\RT WhiteA\RTWhiteAFORC\RT WhiteA 0Hz 7V 1Average Table2.tsv"
    t = 255E-7 
    a = 1E-4 # mask defined area that was used in measurement 
    aReal = 1E-4 # includes effect of undercut during M1 etch
    
### FeFETD1 - FE ###
if device == 2:
    freqdir = r"..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_freqs"
    tempdir = r"..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_temps"
    templkgdir = r"..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_lkg"
    forcFile = r"..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_forc\FeFETD1_die84_MFS+_100_10x10 0Hz 5V 1Average Table1.tsv"
    t = 10E-7 
    a = 1E-4 # mask defined area that was used in measurement 
    aReal = 8.1E3 # includes effect of undercut during M1 etch
###############

#### FeFETD5 - AFE ###
if device == 3:
    freqdir = r"..\tests\testData\FeFETD5\MFS+\die84\FeFETD5_die84_MFS+_60_20x20_freqs"
    tempdir = r"..\tests\testData\FeFETD5\MFS+\die84\FeFETD5_die84_MFS+_60_20x20_temps"
    templkgdir = r"..\tests\testData\FeFETD5\MFS+\die84\FeFETD5_die84_MFS+_60_20x20_lkg"
    forcFile = r"..\tests\testData\FeFETD5\MFS+\die84\FeFETD5_die84_MFS+_60_20x20_FORC\FeFETD5_die68_MFS+_60_20x20_FORC_5V 0Hz 5V 1Average Table2.tsv"
    t = 10E-7
    a = 2.4E-4
    aReal = 2.166E4 # includes effect of undercut during M1 etch
################

landau = lf.LandauFull(thickness=t, area=aReal)
templkgfiles = hd.dir_read(templkgdir)

if tempdir != None:
    tempfiles = hd.dir_read(tempdir)

    if leakageComp:
        tempData = hd.list_read(tempfiles, templkgfiles, plot = False,
                               thickness = t, area = a)
    else:
        tempData = hd.list_read(tempfiles, plot = False,
                                thickness = t, area = a)     
    
    landau.a0 = landau.a0_calc(tempData)

freqfiles = hd.dir_read(freqdir)
freqData = hd.list_read(freqfiles, thickness = t, area = a)

cCompData = freqData[1]

landau.c = landau.c_calc(freqData, plot=1)
compensatedData, landau.pr = landau.c_compensation(cCompData)
hd.hyst_plot([cCompData, compensatedData],
             ['Before', 'After'],
             plot_e=False)
freqCompData = list(map(lambda x:landau.c_compensation(x)[0], freqData))
landau.rho_calc(freqData)



freqDataLkgComp = hd.list_read(freqfiles, templkgfiles, thickness = t, area = a)
cCompDataLkgComp = freqDataLkgComp[1]
hd.hyst_plot([cCompData, cCompDataLkgComp],
             ["With Leakage","Without Leakage"], plot_e=False)

### FORC Calculation


landau_forc = hd.HysteresisData(thickness = t, area = a)
landau_forc.tsv_read(forcFile)
landau_forc.hyst_plot(plot_e=1)
e, er, probs = landau_forc.forc_calc(plot = True)
    
domains = landau.domain_gen(e, er, probs, n=100, plot = False)

elimit = 1.1*max(cCompData.voltage)/t

esweep = np.linspace(-elimit,elimit,num=1000)
esweep = np.append(esweep,esweep[::-1])
res = landau.calc_efe_preisach(esweep, domains, c_add= True, plot=0)

# Plots FORC results vs actual hysteresis measurement ##
fig = plt.figure()
fig.clf()
fig.set_facecolor('white')
ax = fig.add_subplot(111)
ax.plot(cCompData.field*1E-6,cCompData.polarization*1E6,esweep*1E-6,res[0]*1E6)
ax.set_ylabel('Polarization Charge ($\mu{}C/cm^2$)')
ax.set_xlabel('Electric Field (MV/cm)')



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


if tempdir != None:
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
    
    # Following code plots a series of diff temp hystdata files on same plot
    # with leakage current subtraction
    
    if leakageComp:
        hystData = []
        legend = []
        for f in tempfiles:
            data = hd.HysteresisData()
            data.tsv_read(f)
            hystData.append(data)
            legend.append(int(data.temp))
        
        legend = sorted(legend)
        hystData = sorted(hystData, key=lambda data: int(data.temp))
        tempData = sorted(tempData, key=lambda data: int(data.temp))
        
        legend = [str(x)+' K' for x in legend]  
        hd.hyst_plot(tempData, legend)
    
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

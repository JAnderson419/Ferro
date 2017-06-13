#!/usr/bin/env python3
"""
Created on Mon Jun 12 10:45:02 2017

@author: Jackson
"""

import matplotlib.pyplot as plt
import numpy as np
from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')


RTfreqDir = r".\testData\RT WhiteA\RTWhiteAFreq"
RTfreqData = hd.dirRead(RTfreqDir)
RTfreq100hz = r".\testData\RT WhiteA\RTWhiteAFreq\RT WhiteA 100Hz 8V 1Average Table1.tsv"

RT100data = hd.HysteresisData()
RT100data.tsvRead(RTfreq100hz)
RT100data.hystPlot()

RTWhiteFilm = lf.LandauSimple(thickness = 255E-7, area=1E-4)
RTWhiteFilm.c = RTWhiteFilm.cCalc(RTfreqData)
RT100compensated, RTWhiteFilm.pr = RTWhiteFilm.cCompensation(RT100data)

RT100compensated.hystPlot()

forcFile = r".\testData\RT WhiteA\RTWhiteAFORC\RT WhiteA 0Hz 7V 1Average Table7.tsv"
RTWhiteAFORC = hd.HysteresisData(area=1E-4, thickness=255E-7)
RTWhiteAFORC.tsvRead(forcFile)
RTWhiteAFORC.hystPlot(plotE=1)
e, er, probs = RTWhiteAFORC.forcCalc(plot = True)
    
domains = RTWhiteFilm.domainGen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-0.28E6,0.28E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
RTWhiteFilm.calcEfePreisach(esweep, domains, plot=1)
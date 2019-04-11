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
RTfreqFiles = hd.dirRead(RTfreqDir)
RTfreqData = hd.listRead(RTfreqFiles)
RTfreq100hz = r".\testData\RT WhiteA\RTWhiteAFreq\RT WhiteA 100Hz 8V 1Average Table1.tsv"

RT100data = hd.HysteresisData()
RT100data.tsv_read(RTfreq100hz)
RT100data.hyst_plot()

RTWhiteFilm = lf.LandauSimple(thickness = 255E-7, area=1E-4)
RTWhiteFilm.c = RTWhiteFilm.c_calc(RTfreqData)
RT100compensated, RTWhiteFilm.pr = RTWhiteFilm.c_compensation(RT100data)

RT100compensated.hyst_plot()

forcFile = r".\testData\RT WhiteA\RTWhiteAFORC\RT WhiteA 0Hz 7V 1Average Table7.tsv"
RTWhiteAFORC = hd.HysteresisData(area=1E-4, thickness=255E-7)
RTWhiteAFORC.tsv_read(forcFile)
RTWhiteAFORC.hyst_plot(plotE=1)
e, er, probs = RTWhiteAFORC.forc_calc(plot = True)
    
domains = RTWhiteFilm.domain_gen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-0.28E6,0.28E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
RTWhiteFilm.calc_efe_preisach(esweep, domains, plot=1)
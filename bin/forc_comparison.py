#!/usr/bin/env python3
"""
Created on Mon Jun 26 17:47:18 2017

@author: Jackson
"""

from context import LandauFilm as lf
from context import HysteresisData as hd

plt.close('all')

forcFile1 = r"..\ferro\tests\testData\RTWhiteB\RTWhiteB_FORC\RTWhiteB 0Hz 5V 1Average Table1.tsv"
forcFile2 = r"..\ferro\tests\testData\RTWhiteB\RTWhiteB_FORC\RTWhiteB 0Hz 5V 1Average Table3.tsv"
t = 255E-7 
a = 1E-4 # mask defined area that was used in measurement 

forc1 = hd.HysteresisData(thickness = t, area = a)
forc1.tsv_read(forcFile1)
forc1.hyst_plot(plotE=1)
forc1.forc_calc(plot = True)

forc2 = hd.HysteresisData(thickness = t, area = a)
forc2.tsv_read(forcFile2)
forc2.hyst_plot(plotE=1)
forc2.forc_calc(plot = True)
#!/usr/bin/env python3
"""
Created on Mon Jun 26 17:47:18 2017

@author: Jackson
"""
from os.path import join, dirname, realpath
from context import models as lf
from context import data as hd


plt.close('all')
testdatadir = join(dirname(dirname(realpath(__file__))), "tests", "testData")

forcFile1 = join(testdatadir, 'RTWhiteB', 'RTWhiteB_FORC',
                 'RTWhiteB 0Hz 5V 1Average Table1.tsv')
forcFile2 = join(testdatadir, 'RTWhiteB', 'RTWhiteB_FORC',
                 'RTWhiteB 0Hz 5V 1Average Table3.tsv')
t = 255E-7 
a = 1E-4 # mask defined area that was used in measurement 

forc1 = hd.HysteresisData(thickness = t, area = a)
forc1.tsv_read(forcFile1)
forc1.hyst_plot(plot_e=1)
forc1.forc_calc(plot = True)

forc2 = hd.HysteresisData(thickness = t, area = a)
forc2.tsv_read(forcFile2)
forc2.hyst_plot(plot_e=1)
forc2.forc_calc(plot = True)
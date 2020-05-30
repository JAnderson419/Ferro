# ---
# jupyter:
#   jupytext:
#     formats: ipynb,pyscript//py:percent,markdown//md
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown] pycharm={"name": "#%% md\n"}
# # Introduction to Ferro Tutorial
# This notebook is an introduction to the ferro package for ferroelectric data analysis and modeling.
# The data analysis is based on hfo2_mfm.py from the ferro package.
# Analyzed is a Metal-Ferroelectric-Metal capacitor consisting of _
#
# ## Table of Contents
# 0. Introduction
# 1. Data Types, Import, and Visualization
# 2. Leakage Current Compensation
# 3. FORC and Presiach Modeling
# 4. Landau Modeling (WIP)
# 5. Miller Modeling (WIP)

# %% pycharm={"name": "#%%\n"}
# %matplotlib inline


import numpy as np
import matplotlib.pyplot as plt
from os.path import join, dirname, realpath
from ferro import data as hd, models as lf

def isnotebook():
    '''
    Check if in jupyter environment.
    From https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    '''
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

if isnotebook():
    current_folder = globals()['_dh'][0]
else:
    current_folder = dirname(dirname(realpath(__file__)))

DATA_ROOT = join(current_folder, "data")

freqdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_freq")
tempdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_S3_temps")
templkgdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_S3_tempslkg")
forcFile = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_forc",
                "H9 die (9,4) 0Hz 4V 1Average Table1.tsv")

# %% [markdown] pycharm={"name": "#%% md\n"}
#

# %% pycharm={"name": "#%%\n"}
templkgfiles = hd.dir_read(templkgdir)
print(templkgfiles)

# %% pycharm={"name": "#%%\n"}
tempfiles = hd.dir_read(tempdir)
tempData = hd.list_read(tempfiles, templkgfiles)
print(tempData)

# %% pycharm={"name": "#%%\n"}
freqfiles = hd.dir_read(freqdir)
freqData = hd.list_read(freqfiles)
print(freqData)

# %% pycharm={"name": "#%%\n"}
hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)
cCompData = freqData[0]

hfo2.c = hfo2.c_calc(freqData, plot=True)

# %% pycharm={"name": "#%%\n"}
compensatedData, hfo2.pr = hfo2.c_compensation(cCompData)
compensatedData.hyst_plot(plot_e=True)
hfo2.rho_calc(freqData)

# %% pycharm={"name": "#%%\n"}
hfo2.a0 = hfo2.a0_calc(tempData)

# %% pycharm={"name": "#%%\n"}
freqDataLkgComp = hd.list_read(freqfiles, templkgfiles)
cCompDataLkgComp = freqDataLkgComp[0]
hd.hyst_plot([cCompData, cCompDataLkgComp],
             ["With Leakage", "Without Leakage"], plot_e=1)

# %% pycharm={"name": "#%%\n"}
### FORC Calculation

hfo2_forc = hd.HysteresisData(area=6579E-8, thickness=13E-7)
hfo2_forc.tsv_read(forcFile)
hfo2_forc.hyst_plot(plot_e=1)
e, er, probs = hfo2_forc.forc_calc(plot = False)

domains = hfo2.domain_gen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-4.5E6,4.5E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
hfo2.calc_efe_preisach(esweep, domains, plot=1)

# %% pycharm={"name": "#%%\n"}
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

# %% pycharm={"name": "#%%\n"}
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

# %% pycharm={"name": "#%%\n"}
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

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
# # FORC Analysis (Work in Progress)
# <a id='top'></a>
# [Back to Table of Contents](0-Intro.ipynb#top)

# %%
# %matplotlib inline

import numpy as np
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
    root_folder = dirname(globals()['_dh'][0])
else:
    root_folder = dirname(realpath(__file__))

DATA_ROOT = join(root_folder, "tests", "testData")
forcFile = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_forc",
                "H9 die (9,4) 0Hz 4V 1Average Table1.tsv")
freqdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_freq")
freqfiles = hd.dir_read(freqdir)
freqData = hd.list_read(freqfiles, thickness = 13E-7, area=6579E-8)

# %% pycharm={"name": "#%%\n"}
### FORC Calculation

hfo2_forc = hd.HysteresisData(area=6579E-8, thickness=13E-7)
hfo2_forc.tsv_read(forcFile)
hfo2_forc.hyst_plot(plot_e=1)


e, er, probs = hfo2_forc.forc_calc(plot = True)

hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)
hfo2.c = hfo2.c_calc(freqData, plot=1)
compensatedData, hfo2.pr = hfo2.c_compensation(freqData[0])
domains = hfo2.domain_gen(e, er, probs, n=100, plot = False)

esweep = np.linspace(-4.5E6,4.5E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
hfo2.calc_efe_preisach(esweep, domains, plot=1)

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
    root_folder = dirname(globals()['_dh'][0])
else:
    root_folder = dirname(realpath(__file__))

DATA_ROOT = join(root_folder, "tests", "testData")
forcFile = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_forc",
                "H9 die (9,4) 0Hz 4V 1Average Table1.tsv")
freqdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_freq")
freqfiles = hd.dir_read(freqdir)
freqData = hd.list_read(freqfiles, thickness = 13E-7, area=6579E-8)

# %% [markdown]
# We will start by loading the FORC measurement file, which is simply a hysteresis measurement with several cycles
# of increasing voltage magnitude and constant dV/dt.

# %% pycharm={"name": "#%%\n"}
hfo2_forc = hd.HysteresisData(area=6579E-8, thickness=13E-7)
hfo2_forc.tsv_read(forcFile)
hfo2_forc.hyst_plot(plot_e=1)

# %% [markdown] pycharm={"name": "#%% md\n"}
# The hysteresis loop on this measurement looks strange due to leakage current present in the sample,
# but it can be seen from the current data that the loop is being repeated sweeping from 3 V
# down to -3 V over several cycles. This is perhaps better seen in the time domain:

# %% pycharm={"name": "#%%\n"}
hfo2_forc.time_plot()

# %% [markdown]
# The leakage current present in this sample will have minimal impact
# on the final extracted distribution since the switching probability is extracted from the mixed
# derivative $\frac{\delta^2P}{\delta E \delta E_r}$.
# `forc_calc()` can be used on the data to return the experimental e, er, and probability distribution.
# This will be passed to a model object to model the film.

# %% pycharm={"name": "#%%\n"}
e, er, probs = hfo2_forc.forc_calc(plot = True)

# %% [markdown] pycharm={"name": "#%% md\n"}
# It is not time to create a model. Currently, this must be done through a LandauFull model object.
# The modeling interface of ferro is less mature than the data handling and can be expected to undergo future change.
#
# Since data is available for this sample, we will use hysteresis measurements taken at several frequencies
# to extract capacitance from $\frac{i}{\delta V/\delta t}$. Capacitive current compensation will then be performed
# to extract the magnitude of overall remnant polarization. If C or Pr are known, the user can define
# them directly: `hfo2.c = known_c_number`
#

# %% pycharm={"name": "#%%\n"}
hfo2 = lf.LandauFull(thickness = 13E-7, area=6579E-8)
hfo2.c = hfo2.c_calc(freqData, plot=1)
compensatedData, hfo2.pr = hfo2.c_compensation(freqData[0])

# %% [markdown] pycharm={"name": "#%% md\n"}
# With all this in place, it is now possible to generate the domains. For this sample, 100 domains are chosen
# as enough to give an estimate of film performance. This particular measurement has poor Er resolution
# due to a measurement point limit on the ferroelectric tester, so more domains do not improve the fit.
# Stay tuned for future improvements in domain probability interpolation, allowing some of this challenge to be overcome.

# %% pycharm={"name": "#%%\n"}
domains = hfo2.domain_gen(e, er, probs, n=100, plot = True)

# %% [markdown] pycharm={"name": "#%% md\n"}
# Finally, a field can be applied to the model to see the domain response:

# %% pycharm={"name": "#%%\n"}
esweep = np.linspace(-4.5E6,4.5E6,num=1000)
esweep = np.append(esweep,esweep[::-1])
plt.plot(esweep)
polarization, final_domain_state = hfo2.calc_efe_preisach(esweep, domains, plot=1)

# %% [markdown] pycharm={"name": "#%% md\n"}
# The returned polarization array gives polarization values for the film corresponding
# to the applied field, while the final_domain_state is an array of 1 or -1 representing the polarization
# state of the domains. This array can be passed to future `calc_efe_preisach()` 
# function calls to provide an initial value for polarization for minor loop studies
#  (if not specified, domains default to -1 state).

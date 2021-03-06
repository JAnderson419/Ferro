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
# <a id='top'></a>
# # 2 - Data Visualization
#
# [Back to Table of Contents](0-Intro.ipynb#top)
#
# In this section, some of the experimental data visualization tools in ferro will be introduced.
# These allow a researcher to verify that data was imported correctly and visualize trends in experimental data
# before attempting model fits.

# %% pycharm={"name": "#%%\n"}
# %matplotlib inline

import pprint
from os.path import join, dirname, realpath
from ferro import data as hd, aixacct as aix

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
    root_folder = dirname(dirname(dirname(realpath(__file__))))

DATA_ROOT = join(root_folder, "tests", "testData")


tempdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_S3_temps")
templkgdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_S3_tempslkg")
freqdir = join(DATA_ROOT, r"hfo2_MFM", "H9_x9y4_1e4_freq")

templkgfiles = hd.dir_read(templkgdir)
tempfiles = hd.dir_read(tempdir)

pp = pprint.PrettyPrinter(indent=2, width=120, depth=4, compact=True)

# %% [markdown] pycharm={"name": "#%% md\n"}
# <a id='pv'></a>
# ## Hysteresis (P-V) Data Visualization

# %% pycharm={"name": "#%%\n"}
tempdata = hd.list_read(tempfiles, thickness=10E-7)
hd.hyst_plot(tempdata)

# %% [markdown] pycharm={"name": "#%% md\n"}
# Selecting one measurement from the group imported, it is easy to view the associated PV loop:

# %% pycharm={"name": "#%%\n"}
data = tempdata[0]
print(data)


# %% pycharm={"name": "#%%\n"}
data.hyst_plot()

# %% [markdown] pycharm={"name": "#%% md\n"}
# It is just as easy to show the x-axis as field vs voltage:

# %% pycharm={"name": "#%%\n"}
data.hyst_plot(plot_e=True)

# %% [markdown] pycharm={"name": "#%% md\n"}
# In this case, the x-axis looks the same due to the units (MV/cm) and thickness (1E-6) of the sample,
# but the field can be seen more clearly if we assume for a moment that our sample thickness was slightly different:

# %% pycharm={"name": "#%%\n"}
data.thickness = 3E-6 # artificially change thickness for demonstration
data.hyst_plot(plot_e=True)
data.thickness = 1E-6 # change thickness back

# %% [markdown] pycharm={"name": "#%% md\n"}
# To see the raw time series data, it is possible to use the time_plot function, which shows the voltage and current
# as a dashed and solid line respectively vs. measurement time:

# %% pycharm={"name": "#%%\n"}
data.time_plot()

# %% [markdown] pycharm={"name": "#%% md\n"}
# <a id='iv'></a>
# ## Leakage Data Visualization and Hysteresis Leakage Compensation
#
#
# This particular sample has a great deal of leakage current, so it may be useful to attempt leakage current subtraction
# before more detailed analysis. For this particular sample, there are a number of leakage data taken at different
# temperatures:

# %% pycharm={"name": "#%%\n"}
pp.pprint(tempfiles)

# %% [markdown]
# To examine one of these, it can be imported and plotted as a LeakageData object:

# %% pycharm={"name": "#%%\n"}
ldata = hd.LeakageData()
ldata.lcm_read(templkgfiles[0])
print(ldata)
ldata.lcm_plot()

# %% [markdown] pycharm={"name": "#%% md\n"}
# We can also attempt to fit a polynomial to the leakage data using lcm_fit, which uses a fifth order polynomial
# with x and y offset terms. The resulting coefficients can be used to compute the polynomial and pass it to lcm_plot to
# view overlaid with the experimental data

# %% pycharm={"name": "#%%\n"}
ldata.lcm_fit()
print(ldata.lcm_parms)
ldata.lcm_plot()

# %% pycharm={"name": "#%%\n"}
data_comp = data.leakage_compensation(ldata)
data_comp.hyst_plot()

# %% [markdown] pycharm={"name": "#%% md\n"}
# This technique works well for samples with stable leakage processes, but will fail if the sample is 
# degrading between measurements. The compensated and uncompensated data can be compared by calling hyst_plot from
# the ferro.data module itself and passing the two HysteresisData objects:

# %% pycharm={"name": "#%%\n"}
hd.hyst_plot([data, data_comp],['Raw','Compensated'])

# %% [markdown] pycharm={"name": "#%% md\n"}
# This process is not hard for a single measurement, but may become tedious for a large collection of data.
# Rather than performing the fit and compensation manually, the user may instead pass a list of leakage data
# measurements to list_read in addition to the hysteresis data files. Ferro will then look for a temperature
# in the filename of the hysteresis data and attempt to match it with a leakage data measurement of the same temperature.

# %% pycharm={"name": "#%%\n"}
compensated_tempdata = hd.list_read(tempfiles, templkgfiles)
hd.hyst_plot(tempdata)
hd.hyst_plot(compensated_tempdata)

# %% [markdown] pycharm={"name": "#%% md\n"}
# <a id='cv'></a>
# ## Capacitive Current Plotting
#
#
#
# Finally, there also exist functions to plot dV/dt and dP/dV of the hysteresis data.
# These are useful for capacitive current fitting, which will be covered in a future tutorial

# %% pycharm={"name": "#%%\n"}
data_comp.dvdt_plot()

# %% pycharm={"name": "#%%\n"}
data_comp.ncv_plot()


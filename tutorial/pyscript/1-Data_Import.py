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

# %% [markdown]
# <a id='top'></a>
# # 1 - Data Import and Visualization
#
# [Back to Table of Contents](0-Intro.ipynb#top)
#
#
# In its current form, ferro supports simple data import from a TSV file as well as native support for aixACCT
# dynamic hysteresis (inc. FORC) and Leakage (I-V) data. The aixACCT parser also understands PUND and Fatigue data,
# but ferro does not yet include native data types for these measurement formats.
#
# In order to explore the different ways data can be imported to ferro,
# we will start by importing the package and defining our data directories.

# %% pycharm={"name": "#%%\n"}
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

pp = pprint.PrettyPrinter(indent=2, width=120, depth=4, compact=True)
pfile = pprint.PrettyPrinter(indent=2, width=120, depth=4, compact=False)

# %% [markdown] pycharm={"name": "#%% md\n"}
# ## TSV import
# <a id='tsv'></a>
#
# TSV import is the simplest method of loading data into ferro, requiring only a table with time,
# voltage, and current data. When using this method, the data loading function will search for the frequency and
# temperature in the filename in the format `r'\d+Hz'` and `r'\d+C'` or `r'\d+K'`- one or more digits followed by a unit.
# If no temperature is specified, a default value of 300K (27C / 80F) is assigned.
#
# For example, the file **H9 die (9,4) 500Hz 4V 1Average Table17.tsv** when imported would automatically be assigned a
# frequency of 500 Hz and a temperature of 300K:

# %% pycharm={"name": "#%%\n"}
data = hd.HysteresisData()
data.tsv_read(join(freqdir, 'H9 die (9,4) 500Hz 4V 1Average Table17.tsv'))
print(data.freq)
print(data.temp)

# %% [markdown] pycharm={"name": "#%% md\n"}
# To see a quick overview of the measurement without inspecting each property, it is also possible to print the object
# directly:

# %% pycharm={"name": "#%%\n"}
print(data)

# %% [markdown]
# To import a number of individual tsv measurement tables at once, ferro also provides dir_read and list_read functions
# for hysteresisData. `dir_read()` will return a list of files in a directory, which can be passed to `list_read()`, which
# will return a list of hysteresisData objects. If a list of both DHM and IV data is passed, `list_read()` will attempt
# leakage current compensation of the hysteresis data, matching a leakage measurement by temperature.

# %% pycharm={"name": "#%%\n"}
templkgfiles = hd.dir_read(templkgdir)
pp.pprint(templkgfiles)

# %% pycharm={"name": "#%%\n"}
tempfiles = hd.dir_read(tempdir)
tempData = hd.list_read(tempfiles, templkgfiles)
pp.pprint(tempData)


# %% [markdown] pycharm={"name": "#%% md\n"}
# ## aixACCT import
# <a id='aixacct'></a>
#
# Native import for DynamicHysteresisResult amd LeakageResult data from aixACCT .dat ASCII export files is also supported.
# This corresponds to DHM and Leakage measurements (if you are not sure of measurement type, open the .dat file in a
# text editor and look at the first line, it will have a header with measurement type as seen below).

# %% pycharm={"name": "#%%\n"}
aix_file = join(DATA_ROOT, r"hfo2_MFM", 'H9_x9y4_1e4_S3_temps.dat')
num_lines = 5
with open(aix_file, encoding='cp1252') as f:
    head = [next(f) for x in range(num_lines)]
pfile.pprint(head)

# %% [markdown] pycharm={"name": "#%% md\n"}
# The direct read in happens in two parts. The first half, `read_tfdata()`, parses the data text file and loads
# all of the information in the file into a dictionary. This dictionary contains all of the data originally stored
# in the .dat file.

# %% pycharm={"name": "#%%\n"}
data_dict = aix.read_tfdata(aix_file)
pp.pprint(data_dict)

# %% [markdown] pycharm={"name": "#%% md\n"}
# The second half of the read in takes this complete set of data, which may include several IV or PV sweeps, and loads
# each measurement into an appropriate sampleData object using tsv_read():

# %% pycharm={"name": "#%%\n"}
data_list = aix.load_tfdata(data_dict)
pp.pprint(data_list)

# %% [markdown] pycharm={"name": "#%% md\n"}
# ## Radiant Technologies Import
# <a id='rt'></a>
#
# Radiant Technologies Vision software supports data export from the bottom right corner of the plot
# window for a measurement task. A txt file export from these testers looks something like the following:

# %% pycharm={"name": "#%%\n"}
rt_file = join(DATA_ROOT, r"Typical_AB_Data_RT", "TypABdata.Hysteresis.2.txt")
num_lines = 10
with open(rt_file, encoding='cp1252') as f:
    head = [next(f) for x in range(num_lines)]
pfile.pprint(head)

# %% [markdown] pycharm={"name": "#%% md\n"}
# This file contains a single measurement, which can be loaded into an individual HysteresisData object.

# %% pycharm={"name": "#%%\n"}
data = hd.HysteresisData()
data.read_RTHyst(rt_file)
print(data)

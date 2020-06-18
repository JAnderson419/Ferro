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
# 0. [Introduction](0-Intro.ipynb#top)
# 1. [Data Types and Import](1-Data_Import.ipynb#top)
#     1. [TSV Import](1-Data_Import.ipynb#tsv)
#     2. [AixACCT Import](1-Data_Import.ipynb#aixacct)
#     3. [Radiant Technologies Import](1-Data_Import.ipynb#rt)
# 2. [Data Visualization and Leakage Current Compensation](2-Data_Visualization.ipynb#top)
#     1. [HysteresisData PV Plotting](2-Data_Visualization.ipynb#pv)
#     2. [LeakageData IV Plotting](2-Data_Visualization.ipynb#iv)
#     3. [Miscellanious HysteresisData Plotting](2-Data_Visualization.ipynb#cv)
# 3. [FORC and Presiach Modeling (WIP)](3-FORC_Modeling.ipynb#top)
# 4. Landau Modeling (WIP)
# 5. Miller Modeling (WIP)
#

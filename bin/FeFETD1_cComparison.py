#!/usr/bin/env python3
"""
Created on Fri May 26 12:50:08 2017

@author: Jackson
"""

from context import models as lf
from context import data as hd
import matplotlib.pyplot as plt
import numpy as np
from os.path import join, dirname, realpath

plt.close('all')

device = 1
t = 10E-7 # cm

testdatadir = join(dirname(dirname(realpath(__file__))), "tests", "testData")


### FeFET D1 - FE ###
if device == 0:

    freqdirs = [
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_100x100_freq'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_155x155_freq'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_200x200_freq'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_60_20x20_freq'),
    #   r'..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_freqs',
    #   r'..\tests\testData\FeFETD1\MFS+\die116\FeFETD1_die116_MFS+_1_100x100_freq',
    #   r'..\tests\testData\FeFETD1\MFS+\die116\FeFETD1_die116_MFS+_60_20x20_freq'
                ]
    
    lkgdirs = [
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_100x100_leakage'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_155x155_leakage'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_1_200x200_leakage'),
        join(testdatadir, 'FeFETD1', 'MFS+', 'die66',
             'FeFETD1_die66_MFS+_60_20x20_leakage'),
    #   r'..\tests\testData\FeFETD1\MFS+\die84\FeFETD1_die84_MFS+_100_10x10_lkg' ,
    #   r'..\tests\testData\FeFETD1\MFS+\die116\FeFETD1_die116_MFS+_1_100x100_leakage',
    #   r'..\tests\testData\FeFETD1\MFS+\die116\FeFETD1_die116_MFS+_60_20x20_leakage'
                ]
    
    a = np.asarray([1E4,24025,4E4,24000])
    p = np.asarray([400,620,800,4800])
    
    # including undercut
    a = np.asarray([9801,23716,39601,21660])
    p = np.asarray([396,616,796,4560])
    
### FeFET D5 - AFE ###   
elif device == 1:
    freqdirs = [
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_1_100x100_freq'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_1_200x200_freq'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_60_20x20_freq'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_100_10x10_freq'),
                ]
    
    lkgdirs = [
                join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_1_100x100_leakage'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_1_200x200_leakage'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_60_20x20_leakage'),
        join(testdatadir, 'FeFETD5', 'MFS+', 'die82',
             'FeFETD5_die82_MFS+_100_10x10_leakage'),
                ]
    a = np.asarray([1E4,4E4,24000,1E4])
    p = np.asarray([400,800,4800,4000])
    
    # including undercut
    a = np.asarray([9801,39601,21660,8100])
    p = np.asarray([396,796,4560,3600])

c = []
for i,f in enumerate(freqdirs):
    dcfiles = hd.dir_read(f)    
    lkgfile = hd.dir_read(lkgdirs[i])
#    data = hd.list_read(dcfiles, lkgfile)
    data = hd.list_read(dcfiles)
    testfilm = lf.LandauFilm()   
    cde = testfilm.c_calc(data, plot=1)
    c.append(cde)

c = np.asarray(c)


##including undercut
#a1 = .81*a
#p1 = .9*p

fit = np.polyfit(p,c/a,1)
fit_fn = np.poly1d(fit)

print(fit)
    
fig5 = plt.figure()
fig5.set_facecolor('white')
plt.clf()
ax5 = fig5.add_subplot(111)
ax5.plot(a,c/a)
plt.xlabel("Area ($um^2$)")
plt.ylabel('C ($F/um^2$)')
#plt.title('$\\rho{}^-$ as Used for FORC Plot')   

fig4 = plt.figure()
fig4.set_facecolor('white')
plt.clf()
ax4 = fig4.add_subplot(111)
ax4.plot(p,c/a,p,fit_fn(p),'k--')
plt.xlabel("Perimeter ($um$)")
plt.ylabel('C ($F/um^2$)')
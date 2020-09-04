import pytest
from ferro import data as hd
from ferro import aixacct as aix
from ferro import models as md
from os.path import join, dirname, realpath
import csv
import numpy as np
import matplotlib.pyplot as plt
import os

'''
T25 = join(dirname(dirname(realpath(__file__))), "tests", "T_25.csv")
a = []
b = []
with open(T25) as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        a.append(float(row[0]))
        b.append(float(row[1]))
def pulse (t, Vsw, tpw, toffset):
    v = Vsw * np.heaviside(t - toffset, 1) - Vsw * np.heaviside(t - toffset - tpw, 1)
    return v
'''
t = np.linspace(0, 1E-3, 1000) #logspace is probably better here for efficiency
'''
def pulse(t, Vsw, tpw, toffset)
	v = Vsw * np.heaviside(t-toffset, 1) - Vsw * np.heaviside(t-toffset-tpw,1)
	return v

def waveform(Vsw, tpw, toffset):
	for t in tpw:
		p = p + pulse(t, Vsw, tpw, toffset)
	return p
'''

#wave = waveform()
#%%
domain = md.DuChenDomain()
film = md.DuChenFilm(2.8)
film.t = np.linspace(1E-7, 1E-1, 300000)
domain.fit_n_lambda(5E-7, 2E-6)
film.probability_calc(2E-6, 5E-7, domain.n, domain.lamb)
film.prob = film.prob / np.sum(film.prob)
V28 = []
i = 0
while i < 20:
    film.switching_sim(domain.pol_state)
    V28.append(domain.pol_state)
    domain.pol_state = []
    i+= 1

domain1 = md.DuChenDomain()
film1 = md.DuChenFilm(2.6)
film1.t = np.linspace(1E-7, 1E-1, 300000)
domain1.fit_n_lambda(3E-6, 8E-6)
film1.probability_calc(8E-6, 3E-6, domain1.n, domain1.lamb)
film1.prob = film1.prob / np.sum(film1.prob)
V26 = []
i = 0

while i < 20:
    film1.switching_sim(domain1.pol_state)
    V26.append(domain1.pol_state)
    domain1.pol_state = []
    i+= 1

domain2 = md.DuChenDomain()
film2 = md.DuChenFilm(2.4)
film2.t = np.linspace(1E-7, 1E-1, 300000)
domain2.fit_n_lambda(9E-6, 4E-5)
film2.probability_calc(4E-5, 9E-6, domain2.n, domain2.lamb)
film2.prob = film2.prob / np.sum(film2.prob)
V24 = []
i = 0

while i < 20:
    film2.switching_sim(domain2.pol_state)
    V24.append(domain2.pol_state)
    domain2.pol_state = []
    i+= 1
''''''
domain3 = md.DuChenDomain()
film3 = md.DuChenFilm(2.2)
film3.t = np.linspace(1E-7, 1E-1, 300000)
domain3.fit_n_lambda(2E-4, 5E-4)
film3.probability_calc(5E-4, 2E-4, domain3.n, domain3.lamb)
film3.prob = film3.prob / np.sum(film3.prob)
V22 = []
i = 0

while i < 20:
    film3.switching_sim(domain3.pol_state)
    V22.append(domain3.pol_state)
    domain3.pol_state = []
    i+= 1


plt.plot(film.t, film.cdf, label='2.8V')
plt.plot(film.t, film1.cdf, label='2.6V')
plt.plot(film.t, film2.cdf, label='2.4V')
plt.plot(film.t, film3.cdf, label='2.2V')
plt.xscale("log")
plt.legend()
plt.title("Cumulative switching probability distribution")
plt.xlabel("tpw[s]")
plt.ylabel("cdf")
plt.show()
'''
v_array = [V28, V26, V24, V22]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
trace_colors = ['red', 'green', 'blue', 'black']
legends = ['2.8V', '2.6V', '2.4V', '2.2V']
legend_lines =[]
for i, v in enumerate(v_array):
    color = trace_colors[i]
    legend = legends[i]
    for j in range(20):
        if j == 0:
            ax.plot(film.t, v[j], c=color, label=legend)
            continue
        ax.plot(film.t, v[j], c=color)

plt.xscale("log")
plt.title("Single domain switching")
plt.xlabel("tpw[s]")
plt.ylabel("Polarization")
ax.legend()
plt.show()
'''
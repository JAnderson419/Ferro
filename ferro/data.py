#!/usr/bin/env python3
"""
HysteresisData class with functions for reading in and plotting hysteresis data.

TODO: move TF-1000 parsing into python, auto-read thickness & freq from metadata

@author: Jackson Anderson, Rochester Institute of Technology jda4923@rit.edu
"""
from warnings import warn
import re
import copy  # used for creating Ilkg compensated copies of exp data
from os import listdir
from os.path import join, isfile, basename
import matplotlib.pyplot as plt
import numpy as np
import csv
from mpldatacursor import datacursor
from scipy.optimize import curve_fit
from scipy import signal
from scipy.ndimage import filters as flt
from scipy.interpolate import griddata


# matplotlib.rcParams.update({'font.size': 16})


def leakage_func(x, a, b, c, d, e, f, g):
    return (
            a * (x - g) ** 5
            + b * (x - g) ** 4
            + c * (x - g) ** 3
            + d * (x - g) ** 2
            + e * (x - g)
            + f
    )


#    return np.sign(x-d)*a*(np.exp(b*(np.abs(x-d))**c)-1)+np.log(e)


def dir_read(path):
    files = []
    r = re.compile(".*\.tsv$")
    files = [
        join(path, f) for f in listdir(path) if (isfile(join(path, f)) and r.match(f))
    ]
    return files


def list_read(files, leakagefiles=None, plot=False, **kwargs):
    """
    Reads in several hysteresis measurements and creates objects for them.
    
    Parameters
    ----------
    
    files : list
        Paths to tsv data files of hysteresis measurements
    
    leakagefiles : list 
        Paths to leakage data for files.
        If none, leakage compensation is not performed.
        
    plot: bool
        Triggers plotting of leakage data with fit.
        
    kwargs: args
        Arguements to pass to HysteresisData()
    
    Returns
    -------
    
    data_list : list
        HysteresisData objects created from files.
    """
    data_list = []
    for f in files:
        data = HysteresisData(**kwargs)
        data.tsv_read(f)
        if leakagefiles:
            r = re.compile(".*(_| )(" + re.escape(str(data.temp)) + '|' + re.escape(str(int(data.temp))) + ")K.*")
            temp_c = str(data.temp - 273)
            temp_c_int = str(int(data.temp - 273))
            r2 = re.compile(".*(_| )(" + re.escape(temp_c) + '|' + re.escape(temp_c_int) + ")C.*")
            no_match = True
            for j in leakagefiles:
                match = r.match(j)
                match2 = r2.match(j)
                if match or match2:
                    no_match = False
                    ldata = LeakageData()
                    ldata.lcm_read(j)
                    ldata.lcm_fit()
                    if plot:
                        ldata.lcm_plot()
                    data = data.leakage_compensation(ldata)
                else:
                    next
            if no_match:
                raise UserWarning(
                    "No leakage file found to match data. Compensation cannot be performed."
                )
        data_list.append(data)

    return data_list


def hyst_plot(data, legend=None, plot_e=False):
    """
    Plots V vs P, V vs I, and time vs V given hysteresis measurement data.
    
    Parameters
    ----------
    data : list 
        HysteresisData objects to plot.
    legend : list 
        str labels corresponding to data.
    plot_e : bool
        If True plots E instead of P.
    
    Returns
    -------
    n/a
    """
    fig1 = plt.figure()
    fig1.set_facecolor("white")
    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212)

    # creates unique color for each item
    colormap = plt.cm.viridis  # uniform greyscale for printing
    #    colormap = plt.cm.nipy_spectral # diverse color for colorblindness
    ax1.set_prop_cycle("c", [colormap(i) for i in np.linspace(0, 0.6, len(data))])
    ax2.set_prop_cycle("c", [colormap(i) for i in np.linspace(0, 0.6, len(data))])
    lines = []
    for d in data:
        if plot_e:
            line = ax1.plot(1e-6 * d.field, 1e6 * d.polarization)
            lines.append(line[0])
            ax1.set_ylabel("Polarization Charge ($\mu{}C/cm^2$)")

            ax2.plot(1e-6 * d.field, 1e6 * d.current)
            ax2.set_xlabel("Electric Field (MV/cm)")
            ax2.set_ylabel("Current ($\mu{}A$)")

        else:
            line = ax1.plot(d.voltage, 1e6 * d.polarization)
            lines.append(line[0])
            ax1.set_ylabel("Polarization Charge ($\mu{}C/cm^2$)")

            ax2.plot(d.voltage, 1e6 * d.current)
            ax2.set_xlabel("Voltage (V)")
            ax2.set_ylabel("Current ($\mu{}A$)")
    if legend:
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        box2 = ax2.get_position()
        ax2.set_position([box2.x0, box2.y0, box2.width * 0.8, box2.height])
        fig1.legend(lines, legend, loc="center right")


def ncv_plot(data, legend=None, plot_e=False):
    """
    Plots dP/dV (capacitance) of PV data.
    
    Parameters
    ----------
    data : list 
        HysteresisData objects to plot.
    legend : list 
        str labels corresponding to data.
    plot_e : bool
        If True plots E instead of P.
    
    Returns
    -------
    n/a
    """
    fig1 = plt.figure()
    fig1.set_facecolor("white")
    ax1 = fig1.add_subplot(111)

    # creates unique color for each item
    colormap = plt.cm.viridis  # uniform greyscale for printing
    #    colormap = plt.cm.nipy_spectral # diverse color for colorblindness
    ax1.set_prop_cycle("c", [colormap(i) for i in np.linspace(0, 0.6, len(data))])
    lines = []
    for d in data:
        ncv = np.diff(d.polarization) / np.diff(d.voltage)
        ncv_v = 0.5 * (d.voltage[1:] + d.voltage[:-1])
        if plot_e:
            line = ax1.plot(1e-6 * ncv_v / d.thickness, 1e6 * ncv)
            lines.append(line[0])
            ax1.set_xlabel("Electric Field (MV/cm)")
            ax1.set_ylabel("Capacitance ($\mu{}F/cm^2$)")

        else:
            line = ax1.plot(ncv_v, 1e6 * ncv)
            lines.append(line[0])
            ax1.set_xlabel("Voltage (V)")
            ax1.set_ylabel("Capacitance ($\mu{}F/cm^2$)")

    if legend:
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        fig1.legend(lines, legend, loc="center right")


def lcm_plot(data, legend=None):
    """
    Plots leakage current as a function of voltage for a list of LeakageData.
    If parameters have been fit to data, will plot modeled leakage as well.
    
    Parameters
    ----------
    data : list 
        LeakageData objects to plot
    legend : list 
        Str labels corresponding to data
    
    Returns
    -------
    n/a
    """

    fig = plt.figure()
    fig.set_facecolor("white")
    plt.cla()
    ax = fig.add_subplot(111)

    # creates unique color for each item
    #    colormap = plt.cm.viridis # uniform greyscale for printing
    colormap = plt.cm.nipy_spectral  # diverse color for colorblindness
    ax.set_prop_cycle("c", [colormap(i) for i in np.linspace(0, 1, len(data))])

    lines = []
    for d in data:
        line = ax.plot(d.lcm_voltage, 1e6 * d.lcm_current)
        lines.append(line[0])
        if d.lcm_parms != []:
            ax.plot(d.lcm_voltage, 1e6 * leakage_func(d.lcm_voltage, *d.lcm_parms))
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Leakage Current ($\mu{}A$)")

    if legend:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        fig.legend(lines, legend, loc="center right")


class SampleData:
    def __init__(self, thickness=13e-7, area=1e-4, temperature=300):
        """
        Parameters
        ----------
        thickness : float
            Thickness of the sample in cm (used to calculate field)
            
        area: float
            Area of the sample in cm^2. This should match the area that was
            used to calculate polarization charge per unit area for the data.
            
        temperature: float
            Temperature (in Kelvin) at which the measurement was taken.
        
        Returns
        -------
        n/a
        """
        self.sample_name = ""
        self.thickness = thickness  # cm
        self.area = area  # cm^2
        self.temp = temperature  # K


class HysteresisData(SampleData):
    def __init__(self, freq=100.0, **kwargs):
        """
        Inherits SampleData. See that class for info on thickness, area, 
        and temperature
        
        Parameters
        ----------
        thickness : float
            Thickness of the sample in cm (used to calculate field)
            
        area: float
            Area of the sample in cm^2. This should match the area that was
            used to calculate polarization charge per unit area for the data.
            
        temperature: float
            Temperature (in Kelvin) at which the measurement was taken.
        
        Returns
        -------
        n/a
        """
        SampleData.__init__(self, **kwargs)
        self.time = []
        self.voltage = []
        self.current = []
        self.polarization = []
        self.capacitance = []
        self.freq = freq  # Hz

    def __str__(self):
        return f'Hysteresis Data, {len(self.voltage)} points, {min(self.voltage):0.2f} to {max(self.voltage):0.2f} V, ' \
               f'{self.freq} Hz, {self.temp}K, Pmax = {1E6*max(self.polarization):0.2f} uC/cm^2'

    def __eq__(self, other):
        if type(self) == type(other):
            if (self.area == other.area and
                    self.thickness == other.thickness and
                    np.array_equal(self.voltage, other.voltage) and
                    np.array_equal(self.polarization, other.polarization) and
                    np.array_equal(self.time, other.time) and
                    np.array_equal(self.current, other.current)):
                return True
            else:
                return False
        else:
            return False

    @property
    def field(self):
        return self.voltage/self.thickness

    @property
    def dt(self):
        return self.time[1] - self.time[0]


    def tsv_read(self, filename, verbose=False):
        """
        Imports TSV measurement data previously parsed by tfDataTSV_v4.pl. 
        For use with TF-1000 hysteresis measurement data.
        
        Parameters
        ----------
        filename : str
            tsv file (with path) to open and parse. Stores data as 
            the associated HysteresisData object's attributes
        
        Returns
        -------
        n/a
        """
        samplenamematch = re.match(r'^(.*) \d+Hz.*', basename(filename))
        if samplenamematch:
            self.sample_name = samplenamematch.group(1)

        r = re.compile(".*(_| )(\d+)C.*")
        try:
            self.temp = float(r.match(filename).group(2)) + 273
        except AttributeError:
            try:
                r = re.compile(".*(\d+)K.*")
                self.temp = float(r.match(filename).group(2))
            except AttributeError:
                self.temp = 300
                if verbose:
                    print("No temperature specified. Defaulting to 300K")
                next

        r = re.compile(".*(_| )(\d+)Hz.*")
        try:
            self.freq = float(r.match(filename).group(2))
        except AttributeError:
            self.freq = 100
            if verbose:
                print("No frequency specified. Defaulting to 100Hz")
            next

        with open(filename, "r") as data:
            headerLine = data.readline()  # skips first line and saves it for later
            data = csv.reader(data, delimiter="\t")

            for row in data:
                if row:  # ignores blank lines
                    self.time.append(row[0])
                    self.voltage.append(row[1])
                    self.current.append(row[3])
                    self.polarization.append(row[4])

        self.time = np.asfarray(self.time)
        self.voltage = np.asfarray(self.voltage)
        self.current = np.asfarray(self.current)  # A
        self.polarization = 1e-6 * np.asfarray(self.polarization)  # C/cm^2

    def leakage_compensation(self, leakage_data):
        """
        Removes leakage current contribution from hysteresis data using 
        model fit to lcm data.

        Parameters
        ----------
        leakage_data: LeakageData
            Object to use in compensation.
        
        Returns
        -------
        comp_data: HysteresisData
            Deep copy of self with leakage current removed.
        """
        ld = leakage_data

        if ld.lcm_parms == []:
            warn(
                "Please run lcm_fit on the Leakage Data before attempting compensation.",
                RuntimeWarning,
            )
            return self

        comp_data = copy.deepcopy(self)  #

        for j, i in enumerate(comp_data.current):
            ilkg = leakage_func(self.voltage[j], *ld.lcm_parms)
            comp_data.current[j] = i - ilkg

        comp_data.current = comp_data.current - np.mean(comp_data.current)

        testpol = np.zeros(len(comp_data.current))
        for i in range(np.size(testpol)):
            if i == 0:
                next
            else:
                testpol[i] = testpol[i - 1] + comp_data.current[i] * self.dt / self.area

        offset = max(testpol) - (max(testpol) - min(testpol)) / 2

        testpol = testpol - offset
        comp_data.polarization = testpol

        return comp_data

    def bandstop_filter(self, y, freqs=[50, 70], plot=False):
        """
        Experimental
        
        Removes high freq noise from measurement data.
        Filter edge set to 2 times measurement switching frequency.
        Includes some code from scipy.signal.iirfilter example code.
        
        
        Parameters
        ----------
        y : array_like 
            The data to be filtered.
        freqs : list 
            The two freqs defining the edge of the bandstop filter.
        plot : bool
            If True, plots filter response
        
        Returns
        -------
        y : complex ndarray 
            The filtered input data (y). Returned in time domain.
        """
        n = len(y)

        f = np.fft.fftfreq(n, self.dt)  # cycles/sec
        f = 2 * np.pi * self.dt * f  # rad/sample

        b, a = signal.butter(1, freqs / (0.5 / self.dt), btype="bandstop")
        w, h = signal.freqz(b, a, f)

        p = signal.filtfilt(b, a, y)

        if plot:
            fig = plt.figure()
            fig.set_facecolor("white")
            ax = fig.add_subplot(111)
            ax.plot(w / np.pi * (0.5 / self.dt), 20 * np.log10(abs(h)), "o")
            ax.set_xscale("log")
            ax.set_xlabel("Frequency [radians / second]")
            ax.set_ylabel("Amplitude [dB]")
            ax.axis((float(self.freq), 0.5 / self.dt, -100, 10))
            ax.grid(which="both", axis="both")

        return p

    def fft_plot(self, y):
        """
        Takes fft of data and plots in frequency domain.
         
        Parameters
        ----------
        y : np array 
            Data to be plotted
        
        Returns
        -------
        n/a
        """
        n = len(self.time)

        pf = np.fft.fft(y)
        tf = np.linspace(0, 1 / (2 * self.dt), n // 2)

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        ax1.set_title(str(self.freq) + " Hz " + str(self.temp) + " K")
        datacursor(ax1.plot(tf, 2.0 / n * np.abs(pf[0: n // 2])))
        ax1.set_xlabel("frequency")

    #        ax1.set_ylabel('Polarization Charge ($\mu{}C/cm^2$)')

    def hyst_plot(self, plot_e=False):
        """
        Plots V vs P, V vs I, and time vs V given hysteresis measurement data.
        
        Parameters
        ----------
        plot_e : bool
            If True plots E instead of P
        """
        if plot_e:
            fig1 = plt.figure()
            fig1.set_facecolor("white")
            ax1 = fig1.add_subplot(211)
            ax1.set_title(str(self.freq) + " Hz " + str(self.temp) + " K")
            datacursor(ax1.plot(1e-6 * self.field, 1e6 * self.polarization))
            ax1.set_ylabel("Polarization Charge ($\mu{}C/cm^2$)")

            ax2 = fig1.add_subplot(212)
            datacursor(ax2.plot(1e-6 * self.field, 1e6 * self.current))
            ax2.set_xlabel("Electric Field (MV/cm)")
            ax2.set_ylabel("Current ($\mu{}A$)")

        else:
            fig1 = plt.figure()
            fig1.set_facecolor("white")
            ax1 = fig1.add_subplot(211)
            ax1.set_title(str(self.freq) + " Hz " + str(self.temp) + " K")
            datacursor(ax1.plot(self.voltage, 1e6 * self.polarization))
            ax1.set_ylabel("Polarization Charge ($\mu{}C/cm^2$)")

            ax2 = fig1.add_subplot(212)
            datacursor(ax2.plot(self.voltage, 1e6 * self.current))
            ax2.set_xlabel("Voltage (V)")
            ax2.set_ylabel("Current ($\mu{}A$)")

    def time_plot(self):
        """
        Plots forced voltage and measured current vs. time on same plot
        """
        fig2 = plt.figure()
        fig2.set_facecolor("white")
        plt.clf()
        ax4 = fig2.add_subplot(111)
        ax4.plot(self.time, self.voltage, "--")
        ax4.set_xlabel("time (s)")
        ax4.set_ylabel("Voltage (V)")
        ax41 = ax4.twinx()
        ax41.plot(self.time, self.current * 1e6)
        ax41.set_ylabel("Current ($\mu{}A$)")

    def dvdt_plot(self):
        """
        Plots abs(dvdt) of a measurement. 
        Can be used to investigate noise in capacitance extraction.
        """

        dvdt = np.abs(np.diff(self.voltage) / self.dt)
        avg = np.mean(dvdt)

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        ax1.set_title(self.freq + " Hz")
        datacursor(ax1.plot(dvdt, "r."))
        ax1.plot([0, len(dvdt)], [avg, avg], "k--", linewidth=2)
        ax1.set_xlabel("count")
        ax1.set_ylabel("abs(dv/dt) (V/s)")

    def ncv_plot(self, plot_e=False):
        """
        Plots dP/dV (capacitance) of PV data.
        
        Parameters
        ----------
        plot_e : bool
            If True plots E instead of P.
        
        Returns
        -------
        n/a
        """
        fig1 = plt.figure()
        fig1.set_facecolor("white")
        ax1 = fig1.add_subplot(111)

        ncv = np.diff(self.polarization) / np.diff(self.voltage)
        ncv_v = 0.5 * (self.voltage[1:] + self.voltage[:-1])
        if plot_e:
            ax1.plot(1e-6 * ncv_v / self.thickness, 1e6 * ncv)
            ax1.set_xlabel("Electric Field (MV/cm)")
            ax1.set_ylabel("Capacitance ($\mu{}F/cm^2$)")

        else:
            ax1.plot(ncv_v, 1e6 * ncv)
            ax1.set_xlabel("Voltage (V)")
            ax1.set_ylabel("Capacitance ($\mu{}F/cm^2$)")

    def forc_calc(self, plot=False, linear=True, filt_iter=None, filt_dim=[1, 1]):
        """
        Finds minima/maxima in voltage data, cuts down data to only reversal 
        curves, interpolates data onto a linear grid using griddata, and 
        makes contour plot for FORC measurements. 
        
        Note that griddata uses the natgrid library. If you are on windows and
        don't have a visual studio c++ compiler for setup.py, pre-built wheels 
        found at http://www.lfd.uci.edu/~gohlke/pythonlibs/#natgrid.
        (pip install the wheel file directly).
        
        You can also use linear interpolation instead.
        
        Parameters
        ----------
            plot : bool
                Turns plotting of results on if set to True
            linear : bool 
                Switches to linear interpolation
            filt_iter : int
                Number of rolling average filter iterations to apply to the 
                mixed partial derivative of p.
            filt_dim : len 2 sequence of ints
                Size of filter for E & Er axes.
            
        Returns
        ----------
            uniform_e : 1d np array
                E values that correspond to prob
            uniform_er : 1d np array
                Er values that correspond to prob
            prob : 2d numpy array
                Prob distribution of grains for given E,Er
        """

        minima_index = []
        maxima_index = []
        for i, v in enumerate(self.voltage):
            if i == 0:
                continue
            elif i == len(self.voltage) - 2:
                break
            if (v < self.voltage[i - 1]) & (v < self.voltage[i + 1]):
                minima_index.append(i)
            if (v > self.voltage[i - 1]) & (v > self.voltage[i + 1]):
                maxima_index.append(i)

        if len(minima_index) < 2:
            raise ValueError(
                "Could not find two minima for FORC calculation. Aborting."
            )

        vr_forc = []
        v_forc = []
        p_forc = []
        for i, v in enumerate(minima_index):
            vr_forc.extend([self.voltage[v - 1]] * (maxima_index[i + 1] - v))
            v_forc.extend(self.voltage[v: maxima_index[i + 1]])
            p_forc.extend(self.polarization[v: maxima_index[i + 1]])

        e_forc = np.asfarray(v_forc) / (self.thickness)  # V/cm
        er_forc = np.asfarray(vr_forc) / (self.thickness)  # V/cm

        uniform_e = np.linspace(e_forc.min(), e_forc.max(), 200)
        uniform_er = np.unique(np.sort(er_forc))
        uniform_v = np.linspace(
            np.asfarray(v_forc).min(), np.asfarray(v_forc).max(), 200
        )
        uniform_vr = np.unique(np.sort(np.asfarray(vr_forc)))
        xi_grid, yi_grid = np.meshgrid(uniform_e, uniform_er)

        if linear:
            grid = griddata(
                (e_forc, er_forc), p_forc, (xi_grid, yi_grid), method="linear"
            )
        else:
            grid = griddata(
                (e_forc, er_forc), p_forc, (xi_grid, yi_grid), method="nearest"
            )
        grid = np.ma.masked_equal(grid, np.NaN)
        dE, dEr = np.gradient(
            grid, uniform_vr[1] - uniform_vr[0], uniform_v[1] - uniform_v[0]
        )
        dEdE, dEdEr = np.gradient(
            dE, uniform_vr[1] - uniform_vr[0], uniform_v[1] - uniform_v[0]
        )
        if filt_iter != None:
            mask = np.ma.getmask(dEdEr)  # store mask for later
            dEdEr = dEdEr.filled(0)  # fill mask to prevent filter NaN errors
            for i in range(filt_iter):
                dEdEr = flt.uniform_filter(dEdEr, filt_dim)
            dEdEr = np.ma.masked_where(mask, dEdEr)  # reapply mask
        prob = abs(dEdEr) / np.sum(abs(dEdEr))  # normalize for prob dist

        if plot:
            fig4 = plt.figure()
            fig4.set_facecolor("white")
            plt.clf()
            ax4 = fig4.add_subplot(111)
            plt.set_cmap("jet")
            #            forc_plot = ax4.contourf(uniform_v,uniform_vr,1E6*-dEdEr,75)
            forc_plot = ax4.contourf(1e-6 * uniform_e, 1e-6 * uniform_er, prob, 75)
            cbar = plt.colorbar(forc_plot)
            #        cbar.formatter.set_scientific(True)
            cbar.update_ticks()
            cbar.set_label(r"Probability")
            plt.ticklabel_format(style="sci", axis="x", scilimits=(-2, 3))
            plt.ticklabel_format(style="sci", axis="y", scilimits=(-2, 3))
            ax4.set_xlabel("E (MV/cm)")
            ax4.set_ylabel("$E_r$ (MV/cm)")
            ax4.set_title("FORC Plot")
            #            ax4.plot([uniform_vr[0],uniform_vr[-1]],[uniform_vr[0],uniform_vr[-1]],'k-')

            fig5 = plt.figure()
            fig5.set_facecolor("white")
            plt.clf()
            ax5 = fig5.add_subplot(111)
            ax5 = plt.scatter(v_forc, np.asfarray(p_forc) * 1e6, c=vr_forc, alpha=0.25)
            plt.xlabel("$V (V)$")
            plt.ylabel("$P_r (\mu{}C/cm^2)$")
            cb = plt.colorbar(ax5)
            cb.set_label("$V_r$")
            plt.title("$\\rho{}^-$ as Used for FORC Plot")

            fig5 = plt.figure()
            fig5.set_facecolor("white")
            plt.clf()
            ax5 = fig5.add_subplot(111)
            ax5 = plt.scatter(
                e_forc * 1e-6, np.asfarray(p_forc) * 1e6, c=er_forc * 1e-6, alpha=0.25
            )
            plt.xlabel("E (MV/cm)")
            plt.ylabel("$P_r (\mu{}C/cm^2)$")
            cb = plt.colorbar(ax5)
            cb.set_label("$E_r$ (MV/cm)")
            plt.title("$\\rho{}^-$ as Used for FORC Plot")

        return uniform_e, uniform_er, prob


class LeakageData(SampleData):
    def __init__(self, **kwargs):
        SampleData.__init__(self, **kwargs)
        self.lcm_voltage = []
        self.lcm_current = []
        self.lcm_parms = []

    def __eq__(self, other):
        if type(self) == type(other):
            if (self.area == other.area and
                    self.thickness == other.thickness and
                    np.array_equal(self.lcm_voltage, other.lcm_voltage) and
                    np.allclose(self.lcm_current, other.lcm_current) # calculated, so need tolerance for float math
            ):
                return True
            else:
                return False
        else:
            return False

    def lcm_read(self, filename):
        """
        Imports TSV measurement data previously parsed by tfDataTSV_v5.pl
        For use with TF-1000 leakage current measurement data
        
        Parameters
        ----------
        filename : str
            tsv file (with path) to open and parse. Stores data as 
            the associated HysteresisData object's attributes
        
        Returns
        -------
        n/a
        """
        samplenamematch = re.match(r'^(.*) \d*s.*', basename(filename))
        if samplenamematch:
            self.sample_name = samplenamematch.group(1)

        r = re.compile(".*(_| )(\d+)C.*")
        try:
            self.temp = int(r.match(filename).group(2)) + 273
        except AttributeError:
            try:
                r = re.compile(".*(_| )(\d+)K.*")
                self.temp = int(r.match(filename).group(2))
            except AttributeError:
                print("No temperature specified. Defaulting to 300K")
                next

        with open(filename, "r") as data:
            headerLine = data.readline()  # skips first line and saves for later
            data = csv.reader(data, delimiter="\t")

            for row in data:
                if row:  # ignores blank lines
                    self.lcm_voltage.append(row[0])
                    self.lcm_current.append(row[1])

        self.lcm_voltage = np.asfarray(self.lcm_voltage)  # V
        self.lcm_current = self.area * 1e-6 * np.asfarray(self.lcm_current)  # A

    def lcm_fit(
            self,
            func=leakage_func,
            init_guess=np.array([2e-10, 2e-10, 0.8e-6, -1e-6, 1e-6, 0, -1]),
    ):
        """
        Attempts to fit parameters to leakage current, stores in hd object
        
        Parameters
        ----------
        func : function 
            Defines eqn to be used to fit data
        init_guess : np array of appropriate length to match func            
            Provides initial values for curve_fit
        Returns
        -------
        n/a
        """
        # FIXME: curve_fit has trouble converging with some data
        self.lcm_parms, pcov = curve_fit(
            func, self.lcm_voltage, self.lcm_current, p0=init_guess
        )
        print("Fit Parms:", self.lcm_parms)
        print("Std Dev:", np.sqrt(np.diag(pcov)))

    def lcm_plot(self, func=leakage_func):
        """ 
        Plots measured leakage current with fit data.
        """
        fig = plt.figure()
        fig.set_facecolor("white")
        plt.cla()
        ax = fig.add_subplot(111)
        #        datacursor(ax.plot(self.lcm_voltage,np.log(np.abs(self.lcm_current))))
        datacursor(ax.plot(self.lcm_voltage, 1e6 * self.lcm_current, "o"))
        if self.lcm_parms != []:
            #            ax.plot(self.lcm_voltage,np.log(np.abs(leakage_func(self.lcm_voltage,*self.lcm_parms))))
            ax.plot(self.lcm_voltage, 1e6 * func(self.lcm_voltage, *self.lcm_parms))
        ax.set_xlabel("Voltage (V)")
        ax.set_ylabel("Leakage Current ($\mu{}A$)")


def main():
    plt.close("all")


if (
        __name__ == "__main__"
):  # Executes main automatically if this file run directly rather than imported
    main()

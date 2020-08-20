#!/usr/bin/env python3
"""
Created on Mon Apr 10 14:53:15 2017

@author: Jackson Anderson, Rochester Institute of Technology jda4923@rit.edu
"""

# import re
import copy  # used for creating C/Ilkg compensated copies of exp data
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import scipy.stats as rv_discrete
import scipy.stats as ss
import scipy.constants as sc
# from scipy.stats import skew
# from scipy.stats import skewnorm
from scipy.optimize import fsolve, minimize, basinhopping, fmin_slsqp
import numpy as np
from mpldatacursor import datacursor
import csv
import math
import os
import random
# from ferro import data as hd
# from mpl_toolkits.mplot3d import Axes3D


class LandauFilm:
    """
    Base class for Landau Modeling of ferroelectric thin films using alpha
    and beta parameters. 
    
    LandauSimple models behaviour at one temperature (solves for alpha rather
    than alpha0 and Curie Temperature). It also asumes a viscosity coefficient
    of 0.
    
    LandauFull implements rho, Tc, and a0 for more detailed analysis.
    """

    def __init__(self, thickness=13e-7, area=6606e-8, c=0, pr=0):
        self.thickness = thickness  # cm
        self.area = area  # cm^2
        self.c = c  # F
        self.pr = pr

    def c_calc(self, hyst_data, plot=False):
        """
        Calculates non-ferroelectric sample capacitance for the landau film 
        given a list of tf1000 DHM measurement CSV files, 
        with the DHM measurements taken at different frequencies
        
        Calculates non-ferroelectric sample capacitance from the slope of 
        i = C dV/dt using the median abs(current) value at different freq
        (will give non-switching current value)
        
        Parameters
        ----------
        hyst_data : array_like of HysteresisData files.
        plot: Boolean
            Toggles display of matplotlib plot with data fit.

        Returns
        -------
        i_fit[0] : float
            Capacitance in farads
        """

        med_i = np.zeros(len(hyst_data))
        freqs = np.zeros(len(hyst_data))
        dvdt = np.zeros(len(hyst_data))

        for i, d in enumerate(hyst_data):
            freqs[i] = d.freq
            dvdt[i] = np.mean(np.abs(np.diff(d.voltage) / d.dt))
            med_i[i] = np.median(np.abs(d.current))

        i_fit = np.polyfit(dvdt, med_i, 1)
        i_fit_fn = np.poly1d(i_fit)
        er = i_fit[0] * self.thickness / (self.area * sc.epsilon_0 * 1e-2)

        if plot:
            fig2 = plt.figure()
            fig2.set_facecolor("white")
            ax2 = fig2.add_subplot(111)
            ax2.set_title("Capacitance = {:0.3e} F, er= {:.3f}".format(i_fit[0], er))
            datacursor(
                ax2.plot(dvdt, med_i * 1e6, "o", dvdt, 1e6 * i_fit_fn(dvdt), "--k")
            )
            ax2.set_xlabel("dV/dt (V/s)")
            ax2.set_ylabel("Current ($\mu{}A$)")

        return i_fit[0]

    def c_compensation(self, data, plot=False):
        """ 
        Calculates Pr value by subtracting out effect of capacitance in PV curve
        
        Parameters
        ----------
        data: HysteresisData object
        plot: Boolean
            Toggles display of matplotlib plot of original and compensated data.

        Returns
        -------
        comp_data: HysteresisData object
            identical to input but with C*dv/dt current subtracted out
            and polarization recalculated from new current
        pr: float
            remnant polarization value
        """

        # TODO: Test with high leakage current samples

        comp_data = copy.deepcopy(data)  #

        dvdt = np.mean(np.abs(np.diff(data.voltage) / data.dt))
        #        dvdt = (data.voltage[1]-data.voltage[0])/data.dt
        icap = self.c * dvdt

        for j, i in enumerate(data.current):
            if abs(i) >= icap:
                comp_data.current[j] = i - np.sign(i) * icap
            else:
                comp_data.current[j] = 0

        testpol = np.zeros(len(data.current))
        for i in range(np.size(testpol)):
            if i == 0:
                next
            else:
                testpol[i] = (
                    testpol[i - 1]
                    + comp_data.current[i] * comp_data.dt / comp_data.area
                )

        pr = (max(testpol) - min(testpol)) / 2
        offset = max(testpol) - (max(testpol) - min(testpol)) / 2

        testpol = testpol - offset
        comp_data.polarization = testpol

        if plot:
            data.hyst_plot()
            comp_data.hyst_plot()

        return comp_data, pr

    def domain_gen(self, e, er, prob, n=100, plot=False, retParms=False):
        """
        Creates N ferroelectric domains with Ebias and Ec based on the given
        FORC probability distribution.
        
        Gd: HfO2 forms columnar grains roughly 1:1 aspect ratio:
            n = Area/thickness^2
        (See Hoffmann et al., 2016, DOI: 10.1002/adfm.201602869)
          
        Parameters
        ----------
        e: 1d array of uniformly spaced field values from FORC  
        er: 1d array of uniformly spaced reverse field values from FORC            
        prob: 2d array of probabilities (0 to 1) from FORC calculation        
        n: int, number of domains        
        plot: bool, triggers plotting of generated parms        
        retParms: bool, triggers return of array of domain parms
            
        Returns
        -------
        domain_list: list containing domain objects created from generated
                    parameters.

        domains: 2d array of length n by 4 containing domain parameters:
            domains[:,0] = E    
            domains[:,1] = Er  
            domains[:,2] = Ec
            domains[:,3] = Ebias
        """

        prob = np.ma.filled(prob, 0)
        probf = np.ndarray.flatten(prob)
        index = np.zeros(len(probf))
        for i, v in enumerate(index):
            if i == 0:
                next
            else:
                index[i] = index[i - 1] + 1
        #        print (probf)
        domain_list = []
        domains = np.zeros([n, 4])
        for i in range(n):
            j = int(np.random.choice(index, p=probf))
            egrain = e[j % prob.shape[1]]
            ergrain = er[j // prob.shape[1]]
            domains[i, 0] = egrain
            domains[i, 1] = ergrain
            domains[i, 2] = (egrain - ergrain) / 2  # ec
            domains[i, 3] = (egrain + ergrain) / 2  # ebias
            gen_domain = LandauDomain(self, self.area / n, domains[i, 2], domains[i, 3])
            domain_list.append(gen_domain)

        if plot:
            fig1 = plt.figure()
            fig1.set_facecolor("white")
            plt.cla()
            ax1 = fig1.add_subplot(111)
            ax1.plot(1e-6 * domains[:, 0], 1e-6 * domains[:, 1], "o")
            ax1.set_title("Randomly Generated Grain Distribution")
            ax1.set_xlabel("E (MV/cm)")
            ax1.set_ylabel("$E_r$ (MV/cm)")

        if retParms:
            return (domain_list, domains)
        else:
            return domain_list

    def get_ufe(self, pvals, domains):
        """
        Sums potential energy of all ferroelectric domains in a film to get
        the overall energy landscape.
          
        Parameters
        ----------
        pvals: np array of polarization values at which to solve Ufe
        domains: list containing all domain objects in the film
            
        Returns
        -------
            uFE: potential energy density landscape of film
        """

        u = np.zeros(len(pvals))
        for i in domains:
            u = u + i.get_ufe(pvals)

        return u

    def calc_efe_preisach(
        self, esweep, domains, init_state=None, plot=False, c_add=False
    ):
        """
        Models domains as simple hystereons using ec, pr
        
        Parameters
        ----------
        esweep: np array of field values for which to calculate Pr
        domains: list containing all domain objects in the film
        init_state: np array containing initial state of hysterons (-1 or 1)
        plot: bool, triggers plotting of generated hysteresis curve
            
        Returns
        -------
        p: np array, polarization charge values for film (C/cm^2)
        state: np array, final state of hysterons
        """
        if init_state == None:
            state = -np.ones(len(domains))
        else:
            state = init_state

        sweepDir = np.gradient(esweep) #sweep direction

        p = np.zeros(len(esweep)) #esweep = applied electric field
        for j, e in enumerate(esweep):
            for i, d in enumerate(domains): #domains in ferroelectric film
                if sweepDir[j] > 0:
                    if e >= d.ec + d.ebias:
                        state[i] = 1
                elif sweepDir[j] < 0:
                    if e <= -d.ec + d.ebias:
                        state[i] = -1
                # Need to sum actual charge rather than charge density
                p[j] = p[j] + d.pr * d.area * state[i] #total charge at the certain electric field
                #pr is remmenant polarization (charge/area)
                # convert back into charge density
        if c_add:
            p = (p + esweep * self.thickness * self.c) / self.area
        else:
            p = p / self.area

        if plot:
            fig1 = plt.figure()
            fig1.set_facecolor("white")
            plt.cla()
            ax1 = fig1.add_subplot(111)
            ax1.set_title("Presiach Modeled Hysteresis")
            datacursor(ax1.plot(esweep * 1e-6, 1e6 * p))
            ax1.set_xlabel("Electric Field (MV/cm)")
            ax1.set_ylabel("Polarization Charge ($\mu{}C/cm^2$)")

        return p, state

    def u_plot(self, pvals, ufe):
        """
        Plots U vs P for landau film.
        
        Parameters
        ----------
        pvals: 1d np array of polarization charge values
        ufe: 1d np array of energy densities calculated at xVals.
            
        Returns
        -------
        n/a
        """

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        datacursor(ax1.plot(pvals * 1e6, ufe))
        ax1.set_xlabel("Polarization Charge, P (uC/cm^2)")
        ax1.set_ylabel("Energy, U")

    def e_plot(self, pvals, efe, ec=None, ebias=0):
        """
        Plots E vs P for landau film.
        
        Parameters
        ----------
        pvals: 1d np array of polarization charge values
        efe: 1d np array of electric field calculated at xVals.
        ec: float, overlays dotted line at ec+ebias and -ec+ebias on plot
        ebias: float, defines ebias. Used only if ec defined
            
        Returns
        -------
        n/a
        """

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        datacursor(ax1.plot(efe * 1e-6, pvals * 1e6))
        ax1.set_ylabel("Polarization Charge, P (uC/cm^2)")
        ax1.set_xlabel("Electric Field, MV/cm")
        if ec != None:
            e1 = (-ec - ebias) * 1e-6
            e2 = (ec - ebias) * 1e-6
            ax1.plot([e1, e1], [np.min(pvals) * 1e6, np.max(pvals) * 1e6], "r--")
            ax1.plot([e2, e2], [np.min(pvals) * 1e6, np.max(pvals) * 1e6], "r--")


class LandauSimple(LandauFilm):
    """
    Simplified implementation of Landau model. See LandauBase for more info.
    """

    def __init__(self, a=0, **kwargs):
        LandauFilm.__init__(self, **kwargs)
        self.a = a

    def a_calc(self, c=None, t=None):
        if c is None:
            c = self.c
        if t is None:
            t = self.thickness
        a = -1 / (4 * c * t)  # cm/F
        return a


class LandauFull(LandauFilm):
    """
    Full implementation of Landau model. See LandauBase for more info.
    """

    def __init__(self, a0=0, T0=0, rho=0, **kwargs):
        LandauFilm.__init__(self, **kwargs)
        self.a0 = a0
        self.T0 = T0
        self.rho = rho

    def rho_calc(self, hyst_data):
        """
        IN DEVELOPMENT - needs further testing
        
        Calculates a viscosity coefficient for the landau film given a list
        of tf1000 DHM measurement CSV files, with the DHM measurements 
        taken at different frequencies

        Parameters
        ----------
        files : array_like of HysteresisData files.
        
        Returns
        -------
        n/a - not implemented yet 
        """
        # TODO: work on improving rho calculation (noise from C, leakage I)
        c_comp = 0
        pmax = []
        for f in hyst_data:
            pmax.append(np.max(f.polarization))

        pmax = np.asfarray(pmax)
        p = np.min(pmax)

        dpdt = np.zeros(len(hyst_data))
        e = np.zeros(len(hyst_data))

        for i, d in enumerate(hyst_data):
            dt = d.time[1] - d.time[0]
            # print(dvdt[i],"\n", d.file_name)
            dvdt = (d.voltage[1] - d.voltage[0]) / dt

            for j, q in enumerate(d.polarization):
                if j == 0:
                    next
                else:
                    if c_comp == 1:
                        dvdt = (d.voltage[j] - d.voltage[j - 1]) / dt
                        d.current[j] = d.current[j] - self.c * dvdt

                if q < p:
                    next
                elif q >= p:
                    dpdt[i] = (q - d.polarization[j + 1]) / dt
                    f = (p - d.polarization[j - 1]) / (q - d.polarization[j - 1])
                    e[i] = d.field[j - 1] * (1 - f) + f * d.field[j]
                    #                    print (d.polarization[j-1],p,q, '    ',d.field[j-1],e[i],d.field[j])
                    break
        #            d.hyst_plot()
        #            dvdt_filtered = d.bandstop_filter(np.diff(d.voltage)/d.dt)
        #            d.fft_plot(np.diff(d.voltage)/d.dt)
        #            d.fft_plot(dvdt_filtered)

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        datacursor(ax1.plot(dpdt, e))
        ax1.set_xlabel("dP/dt")
        ax1.set_ylabel("dE")

    def a0_calc(self, hystData):
        """
        IN DEVELOPMENT - needs further testing
        
        Calculates an a0 for the landau film given a list
        of tf1000 DHM measurement CSV files, with the DHM measurements 
        taken at different frequencies  
        
        Parameters
        ----------
        files : array_like of HysteresisData files.
        
        Returns
        -------
        n/a - not implemented yet       
        """
        # FIXME: a0 is an order of magnitude too high - need better temp data?
        c_comp = 0
        pmax = []
        for f in hystData:
            pmax.append(np.max(f.polarization))

        pmax = np.asfarray(pmax)
        p = np.min(pmax)

        dpdt = np.zeros(len(hystData))
        e = np.zeros(len(hystData))
        temp = np.zeros(len(hystData))

        for i, d in enumerate(hystData):
            temp[i] = d.temp
            dt = d.time[1] - d.time[0]
            # print(dvdt[i],"\n", d.file_name)
            dvdt = (d.voltage[1] - d.voltage[0]) / dt

            for j, q in enumerate(d.polarization):
                if j == 0:
                    next
                else:
                    if c_comp == 1:
                        dvdt = (d.voltage[j] - d.voltage[j - 1]) / dt
                        d.current[j] = d.current[j] - self.c * dvdt

                if q < p:
                    next
                elif q >= p:
                    dpdt[i] = (q - d.polarization[j + 1]) / dt
                    f = (p - d.polarization[j - 1]) / (q - d.polarization[j - 1])
                    e[i] = d.field[j - 1] * (1 - f) + f * d.field[j]  # linear interp
                    #                    print (d.polarization[j-1],p,q, '    ',d.field[j-1],e[i],d.field[j])
                    break
        #            d.hyst_plot()
        #            dvdt_filtered = d.bandstop_filter(np.diff(d.voltage)/d.dt)
        #            d.fft_plot(np.diff(d.voltage)/d.dt)
        #            d.fft_plot(dvdt_filtered)

        a0fit = np.polyfit(temp, e, 1)
        a0fit_fn = np.poly1d(a0fit)
        a0 = a0fit[0] / (2 * p)  # cm/(F*K)
        Tc = 1 / (4 * self.c * a0 * self.thickness) + 300  # 300K temp at which C calced
        #        print (a0,Tc)

        fig1 = plt.figure()
        fig1.set_facecolor("white")
        plt.cla()
        ax1 = fig1.add_subplot(111)
        datacursor(ax1.plot(temp, e, "o", temp, a0fit_fn(temp)))
        ax1.set_title("a0 = {:0.3e}".format(a0))
        ax1.set_xlabel("T (C)")
        ax1.set_ylabel("E (V/cm)")


class LandauDomain:
    """
    IN DEVELOPMENT - includes different experimental solving methods for 
        landau parameters
    
    Represents individual ferroelectric domains in the film.
    Used to calculate domain-specific beta parameter for multidomain simulation.
    aTerm is alpha or alpha*(T-Tc) coefficient, depending on film model used.
    """

    def __init__(self, landau, area, ec, ebias, a_term=0, b=0, g=0):
        self.pr = landau.pr
        self.t = landau.thickness
        self.c = landau.c
        self.ec = ec
        self.ebias = ebias
        self.a = a_term
        self.b = b
        self.g = g
        self.area = area

    def eqns(self, p):
        a, b, g = p
        #        print (a,b,g,self.c,self.t,self.pr,self.ec,self.ebias)
        eq1 = (
            self.ec
            - self.ebias
            + 1 * a * self.pr
            + 1 * b * self.pr ** 3
            + 1 * g * self.pr ** 5
        )
        eq2 = (
            2 * a * self.pr
            + 4 * b * self.pr ** 3
            + 6 * g * self.pr ** 5
            + self.ebias
            + self.ec
        )
        eq3 = (2 * a + 12 * b * self.pr ** 2 + 30 * g * self.pr ** 4) - 1 / (
            self.c * self.t
        ) / self.area
        #        print (eq1, eq2, eq3)
        eq4 = 2 * a + 12 * b * self.pr ** 2 + 30 * g * self.pr ** 4

    #        if (a<0 and b>0 and g>0):
    #            return(eq1, eq2, eq3)
    #        else:
    #            return(1E4+a,1E4+a,1E4+a)
    #        if(eq4 > 0 and a<0):
    #            return(eq1, eq2, eq3)
    #        else:
    #            return(1E8,1E8,1E8)
    #        return(eq1, eq2, eq3)

    def eqns1(self, p):
        a, b, g = p
        #        print (a,b,g,self.pr,self.ec,self.ebias)
        eq1 = (
            self.ec
            - self.ebias
            + 1 * a * self.pr
            + 1 * b * self.pr ** 3
            + 1 * g * self.pr ** 5
        ) ** 2
        eq2 = (
            2 * a * self.pr
            + 4 * b * self.pr ** 3
            + 6 * g * self.pr ** 5
            + self.ebias
            + self.ec
        ) ** 2
        eq3 = (
            (2 * a + 12 * b * self.pr ** 2 + 30 * g * self.pr ** 4)
            - 1 / (self.c * self.t) / self.area
        ) ** 2
        #        print (eq1, eq2, eq3)

        return (eq1, eq2, eq3)

    def con(self, p):
        """
        Defines parameter constraints.
        """
        a, b, g = p
        eq4 = 2 * a + 12 * b * self.pr ** 2 + 30 * g * self.pr ** 4
        return eq4

    def parm_calc(self):
        guess = np.asarray(
            (-2 / self.pr ** 2, 1 / (2 * self.pr ** 4), 1 / self.pr ** 6)
        )
        #        self.a,self.b,self.g = fsolve(self.eqns, guess)

        #        print (guess)
        parms = minimize(
            self.eqns1,
            x0=np.asarray(
                (-2 / self.pr ** 2, 1 / (2 * self.pr ** 4), 1 / self.pr ** 6)
            ),
            bounds=((None, 0), (0, None), (0, None)),
            #                         constraints = ({'type':'ineq', 'fun':self.con}),
            method="COBYLA",
        )
        #        kwargs = {"method": "COBYLA",
        #                  "bounds" : ((None,0),(0,None),(0,None)),
        #                  "constraints" : ({'type':'ineq', 'fun':self.con})}
        #        parms = basinhopping(self.eqns1,
        #                             x0 = guess,
        #                             niter = 100,
        #                             minimizer_kwargs = kwargs,
        #                             stepsize = self.pr,
        #                             T = 1/self.pr**2)
        self.a, self.b, self.g = parms.x

    def parm_fit(self, plot=False):
        guess = np.asarray(
            (-2 / self.pr ** 2, 1 / (2 * self.pr ** 4), 1 / self.pr ** 6)
        )
        a_guess, b_guess, g_guess = guess
        guessRes = 26

        avals = np.linspace(a_guess / 1e3, a_guess * 1e3, guessRes)
        bvals = np.linspace(b_guess / 1e3, b_guess * 1e3, guessRes)
        gvals = np.linspace(g_guess / 1e4, g_guess * 1e4, guessRes)

        err = np.empty([guessRes, guessRes, guessRes])
        for i, a in enumerate(avals):
            for j, b in enumerate(bvals):
                vals = np.asfarray(self.eqns1((a, b, gvals)))
                err[i, j, :] = np.sum(vals ** 2, 0)

        if plot:
            x = np.empty(guessRes ** 3)
            y = np.empty(guessRes ** 3)
            z = np.empty(guessRes ** 3)
            d = np.empty(guessRes ** 3)
            for i, r in enumerate(err.flatten()):
                xcord = i % guessRes
                ycord = (i // guessRes) % guessRes
                zcord = i // guessRes ** 2

                x[i] = avals[xcord]
                y[i] = bvals[ycord]
                z[i] = gvals[zcord]
                d[i] = r

            fig1 = plt.figure()
            #        colormap = plt.cm.viridis # uniform greyscale for printing
            colormap = plt.cm.nipy_spectral  # diverse color for colorblindness
            plt.clf()
            ax1 = fig1.add_subplot(111, projection="3d")
            #            ax1.set_xscale('log')
            #            ax1.set_yscale('log')
            #            ax1.set_zscale('log')
            p = ax1.scatter(
                np.abs(x),
                y,
                z,
                c=d,
                alpha=0.5,
                s=15,
                lw=0,
                cmap=colormap,
                norm=Normalize(),
            )
            ax1.set_xlabel(r"-$\alpha{}$ Guess")
            ax1.set_ylabel(r"$\beta{}$ Guess")
            ax1.set_zlabel(r"$\gamma{}$ Guess")
            plt.colorbar(p, ax=ax1)

    def get_ufe(self, pvals):
        """
        Parameters
        ----------
        pvals: np array 
            p values for which to solve ufe
        
        Returns
        -------
        np array
            ufe values
        """
        return (
            self.a * pvals ** 2
            + self.b * pvals ** 4
            + self.g * pvals ** 6
            - self.ebias * pvals
        )

    def get_efe(self, pvals):
        """
        Parameters
        ----------
        pvals: np array 
            p values for which to solve efe
        
        Returns
        -------
        np array 
            efe values
        """
        return (
            2 * self.a * pvals
            + 4 * self.b * pvals ** 3
            + 6 * self.g * pvals ** 5
            - self.ebias
        )


class DuChenDomain():
    def __init__(self, t0=None, alpha=None):
        self.alpha = alpha
        self.tao = 0
        self.phi = 0
        self.t0 = t0
        self.n = 0
        self.lamb = 0
        self.pol_state = []

    def fit_t0_alpha(self, filename, vsw):
        """
        *TESTING REQUIRED*
        Find t0 and alpha based on temperature data
        Parameters
        ----------
        filename: 
            The temperature data with average switching time
        vsw:
            Signal amplitude
        Returns
        -------
        """ 
        with open(filename) as f:
            csv_reader = csv.reader(f, delimiter=',')
            vsw_pow_neg2 = []
            tsw_pow_neg1 = []
            for row in csv_reader:
                vsw_pow_neg2.append(float(row[0]))
                tsw_pow_neg1.append(float(row[1]))
            self.t0 = 1 / (np.polyfit(tsw_pow_neg1, vsw_pow_neg2, 1))[1]
            self.alpha = pow(vsw, 2) * 1.38e-23 * np.log(self.tao / self.t0)

    def fit_n_lambda(self, stdev, mean):
        """
        Find lambda and n based on average and standard deviation switching time
        Parameters
        ----------
        stdev:
            Standard Deviation of switching time
        mean:
            Average switching time

        Returns
        -------
        """
        self.lamb = mean / pow(stdev, 2)
        self.n = mean * self.lamb
        self.tao = 1/self.lamb

class DuChenFilm():
    def __init__(self, vsw):
        #self.temperature = temp
        #self.area = area
        self.polarization = []
        self.vsw = vsw
        self.voltage = []
        self.t = []
        self.tpw = []
        self.domains = []
        self.tpw_avg = 0
        self.tpw_stdev = 0
        self.cdf = []
        self.pdf =[]
        self.prob = []

    def probability_calc(self, mean, stdev, n, lamb):
        """
        Getting cumulative, probability, and gamma distribution based on given parameters
        Parameters
        ----------
        mean:
            average switching time
        stdev:
            standard deviation of switching time
        n:
            number of nucleation (n)
        lamb:
            rate of nucleation (lambda)

        Returns
        -------

        """
        self.prob = lamb * np.exp(-lamb * self.t) * (lamb * self.t) ** (n - 1) / math.factorial(int(n - 1))
        self.cdf = ss.norm.cdf(self.t, mean, stdev)
        self.pdf = ss.norm.pdf(self.t, mean, stdev)

    def switching_sim(self, state):
        """
        Simulate domain switching based on gamma distribution stored in film class
        Parameters
        ----------
        state:
            domain state
        Returns
        -------
        """
        initial_state = -1
        tsw = np.random.choice(self.t, p=self.prob)
        for i, time in enumerate(self.t):
            if self.t[i] < tsw:
                state.append(initial_state)
                continue
            state.append(1)
    #Multidomain films: bias to film
    #Each domain in film have different tsw and vsw relationship

    def read_exp_data(self, expdir):
        """
        *TESTING REQUIRED*
        Read experimental data (voltage, time, and polarization changes)
        Parameters
        ----------
        expdir: experiment data directory

        Returns
        -------
        """
        entries = os.listdir(expdir)
        for entry in entries:
            csv_reader = csv.reader(f, delimiter=',')
            with open(entry) as f:
                i = 0
                time = []
                voltage = []
                polarization = []
                for row in csv_reader:
                    if i != 0:
                        time.append(float(row[0]))
                        voltage.append(float(row[1]))
                        polarization.append(float(row[2]))
                    i += 1
                self.t.append(time)
                self.voltage.append(voltage)
                self.polarization.append(polarization)
                time.clear()
                voltage.clear()
                polarization.clear()

    def extract_tpw(self):
        """
        *TESTING REQUIRED*
        Get switching time based on ideal signal and polarization (need more works)
        Returns
        -------

        """
        for i, element in enumerate(self.polarization):
            dx = self.t[i][1] - self.t[i][0]
            dy = np.diff(element) / dx
            tsw = self.t[i][dy.index(max(dy))]
            index = self.t[i](tsw)
            k = index
            while self.voltage[i][k] == self.voltage[i][index]:
                k -= 1
            self.tpw.append(self.t[i][index] - self.t[i][k+1])
        self.tpw_avg = sum(self.tpw) / len(self.tpw)
        self.tpw_stdev = np.std(self.tpw)

# 10uC/cm^2 = polarization/area)
# Polarization =
# Initial polariza
def main():
    plt.close("all")


#    FeFETD2 = HysteresisData()
#    FeFETD2.tsv_read(r'D:\Google Drive\Ferroelectric Research\FE_20152016\SeniorDesign\Testing\FETester\NaMLab_FeFETD2_die24_100x100_450V_data.csv')
#    FeFETD2.hyst_plot()

#    FeFETD2 = HysteresisData()
#    FeFETD2.tsv_read(r'D:\Google Drive\Ferroelectric Research\FE_20142015\TF-1000\Measurements\RadiantFEcaps_MIM_DiscreteCap\tool_tests\10uF discrete cap 100Hz 5V.tsv')
#    FeFETD2.hyst_plot()

if (
    __name__ == "__main__"
):  # Executes main automatically if this file run directly rather than imported
    main()

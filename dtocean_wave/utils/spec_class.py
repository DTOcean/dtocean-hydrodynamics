# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2017-2018 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
This module contains the classes used to generate the wave spectrums associated with the different sea states

.. module:: spec_class
   :platform: Windows
   :synopsis: Class to generate wave models

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt
from WatWaves import *
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D

style.use("ggplot")


class wave_spec():
    """
    wave_spec: the class collects methods to assess the power spectrum density for the given sea states

    Args:
        f (numpy.ndarray) [Hz]: wave frequencies
        fp (float) [Hz]: peak wave frequency
        Hs (float): significant wave height

    Optional args:
        t (numpy.ndarray) [rad]: wave directions
        t_mean (float) [rad]: mean wave direction
        s (float): direction spreading parameter
        st (str): spectral shape name "Regular",
                                        "Bretschneider_Mitsuyasu",
                                        "Modified_Bretschneider_Mitsuyasu",
                                        "Pierson_Moskowitz",
                                        "Jonswap",
                                        "pscSwell"
        gamma (float): peak enanchement factor
        d (float): water depth
        COFact_Low (float): cut off frequency factor in the lower range. Relative to the peak frequency
        COFact_Top (float): cut off frequency factor in the upper range. Relative to the peak frequency
        COFreq_Low (float): cut off frequency in the lower range
        COFreq_Top (float): cut off frequency in the upper range

    Attributes:
        spec_type (str): spectral shape name "Regular",
                                        "Bretschneider_Mitsuyasu",
                                        "Modified_Bretschneider_Mitsuyasu",
                                        "Pierson_Moskowitz",
                                        "Jonswap",
                                        "pscSwell"
        f (numpy.ndarray) [Hz]: wave frequencies
        fp (float) [Hz]: peak wave frequency
        Hs (float): significant wave height
        gamma (float): peak enanchement factor
        t (numpy.ndarray) [rad]: wave directions
        t_mean (float) [rad]: mean wave direction
        s (float): direction spreading parameter
        d (float): water depth
        CutOffFactorLowBound (float): cut off frequency factor in the lower range. Relative to the peak frequency
        CutOffFactorTopBound (float): cut off frequency factor in the upper range. Relative to the peak frequency
        CutOffFreqLowBound (float): cut off frequency in the lower range
        CutOffFreqTopBound (float): cut off frequency in the upper range
        g (float): acceleration of gravity
        nfrqs (int): number of frequencies
        specs (list) [m^2/Hz/rad]: list of power spectral density distributions
        spec_type_list (list): list of spectral types
        df (Hz): frequency discretisation step
        dth (rad): angular discretisation step
    """
    def __init__(self, f, fp, Hs,
                 t=np.array([0.]), 
                 t_mean=0.,
                 s=0.,
                 st="",
                 gamma=3.3,
                 d=100.,
                 COFact_Low=0.0033,
                 COFact_Top=100.33,
                 COFreq_Low=-1.,
                 COFreq_Top=1e10):
                     
        self.spec_type = st
        self.f = f
        self.fp = fp
        self.Hs = Hs
        self.gamma = gamma
        self.d= d
        self.CutOffFactorLowBound = COFact_Low
        self.CutOffFactorTopBound = COFact_Top
        self.CutOffFreqLowBound = COFreq_Low
        self.CutOffFreqTopBound = COFreq_Top
        self.g = 9.82
        self.nfrqs = len(f)
        self.specs = []
        self.specs_type_list = []
        self.df = 1.
        self.dth = 1.
        
        if np.shape(t)==(1,):
            self.t = np.array([t_mean],'f')
        elif np.shape(t)==():
            errStr('The format of the direction vector is incorrect')
            raise IOError(errStr)
        else:            
            self.t = t
        
        if self.f.shape[0] > 1:
            self.df = np.abs(f[0]-f[1])
            
        if self.t.shape[0] > 1:
            self.dth = np.abs(t[0]-t[1])
            
        self.t_mean = t_mean
        self.s = s
        
        
        if not st=="":
            add_spectrum(self)
            
    def show(self, index = []):
        """
        show: visualises the spectral shapes in function of frequencies and direcitons

        Args:
            index: index of the spectral shape to be visualised
        """
        if index == []:
            index = range(len(self.specs))
        
        fig = plt.figure()

        a = fig.gca(projection='3d')
        a.set_xlabel('frequency, [Hz]')
        a.set_zlabel('psd, [-]')
        a.set_ylabel('direction, [degree]')
        
        specNum = 0
        selectedSpectrum = [self.specs[el] for el in index ]
        for indSp, specs in enumerate(selectedSpectrum):
            specNum += 1
            X,Y = np.meshgrid(specs[0],specs[1])
            a.plot_wireframe(X.T,Y.T*180./np.pi,specs[2],label='Spectrum #{}\nType: {}'.format(specNum,self.specs_type_list[indSp]))

        plt.legend()
        plt.show()
        
    def add_spectrum(self):
        """
        add_spectrum: calls required function to generate the spectrum and adds the current spectrum to the spectrum list

        """
        f = 0.
        S = 0.
        if self.spec_type == "Regular":
            f,S = self.Regular()
            self.s = 0.
        elif self.spec_type == "Bretschneider_Mitsuyasu":
            f,S = self.Bretschneider_Mitsuyasu()
        elif self.spec_type == "Modified_Bretschneider_Mitsuyasu":
            f,S = self.Modified_Bretschneider_Mitsuyasu()
        elif self.spec_type == "Pierson_Moskowitz":
            f,S = self.Pierson_Moskowitz()
        elif self.spec_type == "Jonswap":
            f,S = self.Jonswap()
        elif self.spec_type == "pscSwell":
            f,S = self.pscSwell()
        else:
            errStr = ("The string input is not valid "
                          "The list of input is given in the documentation.")
            raise ValueError(errStr)  
        
        self.specs.append([f,self.t,self.Directional(S)])
        self.specs_type_list.append(self.spec_type)

    def rm_spectrum(self,index=0):
        """
        rm_spectrum: remove the specified spectrum from the spectrum list

        Optional args:
            index (int/str): index of the spectrum to be removed from the list. If "all", the spectrum list will be refreshed
        """
        if index=='all':
             self.specs = []
             self.specs_type_list = []
        else:
            self.specs.pop(index)
            self.specs_type_list.pop(index)
            
    def Regular(self):
        """
        Regular: defines a wave spectrum with dirac delta distribution

        Returns:
            (list): frequencies and spectral distribution
        """
        module_logger.info('WARNING!!! The sqrt(2) factor between the') 
        module_logger.info('significant wave height and the wave height of a regular wave')
        module_logger.info('is not considered!!')
        T = 1/self.fp
        H = self.Hs
        f = self.f
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        
        NFFT_ny = len(f)
        fny = f[-1]
        fs = fny*2
        NFFT = NFFT_ny*2
        dt = 1/fs
        t = np.linspace(0,NFFT,NFFT,endpoint=False)*dt
        eta = H/2.*np.cos(2.*np.pi/T*t)
        ETA = np.fft.fft(eta)
        S_ss = 2.*(f[1]-f[0])*np.abs(ETA[0:NFFT_ny])
        spec_v = np.zeros(NFFT_ny)
        index = (f > np.max((CutOffFreqLowBound,
                             self.fp * CutOffFactorLowBound))) * \
                (f < np.min((CutOffFreqTopBound,
                             self.fp * CutOffFactorTopBound)))
        
        for ifr in range(NFFT_ny):
            if index[ifr]:
                spec_v[ifr] = -S_ss[ifr]
        
        return [f,spec_v]
                
    def Bretschneider_Mitsuyasu(self):
        """
        Bretschneider_Mitsuyasu: defines a wave spectrum with spectral shape defined by the Bretschneider_Mitsuyasu model

        Returns:
            (list): frequencies and spectral distribution
        """

        spec_type = self.spec_type
        f = self.f
        fp = self.fp
        Hs = self.Hs
        gamma = self.gamma
        d = self.d
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        g = self.g
        nfrqs = self.nfrqs
        spec_v = np.zeros((nfrqs))
        
        index = (f > np.max((CutOffFreqLowBound,fp*CutOffFactorLowBound))) * (f<np.min((CutOffFreqTopBound,fp*CutOffFactorTopBound)))
        for ifr in range(nfrqs):     
            if index[ifr]:
                fr = f[ifr]
                spec_ = 1.03*(fp/fr)**4
                if spec_ < 50:
                    spec_v[ifr] = 0.257*Hs**2*fp**4/fr**5 * np.exp(-spec_)
        return [f,spec_v]
        
    def Modified_Bretschneider_Mitsuyasu(self):
        """
        Modified_Bretschneider_Mitsuyasu: defines a wave spectrum with spectral shape defined by the Modified_Bretschneider_Mitsuyasu model

        Returns:
            (list): frequencies and spectral distribution
        """
        spec_type = self.spec_type
        f = self.f
        fp = self.fp
        Hs = self.Hs
        gamma = self.gamma
        d = self.d
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        g = self.g
        nfrqs = self.nfrqs
        spec_v = np.zeros((nfrqs))
        
        index = (f > np.max((CutOffFreqLowBound,fp*CutOffFactorLowBound))) * (f<np.min((CutOffFreqTopBound,fp*CutOffFactorTopBound)))
        for ifr in range(nfrqs):
                 if index[ifr]:
                     fr = f[ifr]
                     spec_ = 0.75*(fp/fr)**4
                     if spec_ < 50:
                         spec_v[ifr] = 0.205*Hs**2*fp**4/fr**5 * np.exp(-spec_)
        return [f,spec_v]
   
    def Pierson_Moskowitz(self):
        """
        Pierson_Moskowitz: defines a wave spectrum with spectral shape defined by the Pierson_Moskowitz model

        Returns:
            (list): frequencies and spectral distribution
        """
        spec_type = self.spec_type
        f = self.f
        fp = self.fp
        Hs = self.Hs
        gamma = self.gamma
        d = self.d
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        g = self.g
        nfrqs = self.nfrqs
        spec_v = np.zeros((nfrqs))
        
        index = (f > np.max((CutOffFreqLowBound,fp*CutOffFactorLowBound))) * (f<np.min((CutOffFreqTopBound,fp*CutOffFactorTopBound)))
        for ifr in range(nfrqs):
           if index[ifr]:
               fr = f[ifr]
               a = 1.25*(fp/fr)**4
               if a<50:
                   spec_v[ifr] = 0.25 *Hs**2*a*np.exp(-a)/fr
               else:
                   spec_v[ifr] =0.

        return [f,spec_v]

    def Jonswap(self):
        """
        Jonswap: defines a wave spectrum with spectral shape defined by the Jonswap model

        Returns:
            (list): frequencies and spectral distribution
        """
        spec_type = self.spec_type
        f = self.f
        fp = self.fp
        Hs = self.Hs
        gamma = self.gamma
        d = self.d
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        g = self.g
        nfrqs = self.nfrqs
        spec_v = np.zeros((nfrqs))
        
        index = (f > np.max((CutOffFreqLowBound,fp*CutOffFactorLowBound))) * (f<np.min((CutOffFreqTopBound,fp*CutOffFactorTopBound)))
        for ifr in range(nfrqs):
           if index[ifr]:
               fr = f[ifr]
               BetaJ = 0.0624/(0.230+0.0336*gamma-0.185*1./(1.9+gamma))
               Tp = 1/fp
               if fr<=fp:
                   Sigma = 0.07
               else:
                   Sigma = 0.09
               spec_v[ifr] = BetaJ * Hs * Hs * Tp**(-4) * fr**(-5) * \
                   np.exp(-1.25 * (Tp*fr)**(-4)) * \
                   gamma**np.exp(-(Tp * fr - 1)**2 / (2. * Sigma**2))
        return [f,spec_v]

    def pscSwell(self):
        """
        pscSwell: defines a wave spectrum with spectral shape defined by the pscSwell model

        Returns:
            (list): frequencies and spectral distribution
        """
        spec_type = self.spec_type
        f = self.f
        fp = self.fp
        Hs = self.Hs
        gamma = self.gamma
        d = self.d
        CutOffFactorLowBound = self.CutOffFactorLowBound
        CutOffFactorTopBound = self.CutOffFactorTopBound
        CutOffFreqLowBound = self.CutOffFreqLowBound
        CutOffFreqTopBound = self.CutOffFreqTopBound
        g = self.g
        nfrqs = self.nfrqs
        spec_v = np.zeros((nfrqs))
        
        index = (f > np.max((CutOffFreqLowBound,fp*CutOffFactorLowBound))) * (f<np.min((CutOffFreqTopBound,fp*CutOffFactorTopBound)))
        for ifr in range(nfrqs):
            if index[ifr]:
                fr = f[ifr]
                spec_ = 1.2/fp/fr*np.sqrt(np.sqrt(fp/fr))
                if spec_>50:
                    spec_ = 0.
                else:
                    spec_ = np.exp(-spec_)
                
                if spec_>0:
                    spec_ = 6. / (2. * np.pi * 16.) * np.sqrt(Hs) * fp * np.sqrt(np.sqrt(np.sqrt(fp / fr))) * spec_
                
                spec_v[ifr] = spec_
                
        return [f,spec_v]

    def Directional(self, S):
        """
        Directional: distribute the spectral energy along the given direction,
        given the spectral parameter

        Args:
            S (numpy.ndarray): power spectral density distribution

        Returns:
            S_d (numpy.ndarray):
                power spectral density distributed along the different
                directions
        """
            
        def LnGamma(X):
            xx = X - 1.0
            tmp = xx + 5.5
            tmp = (xx + 0.5) * np.log(tmp) - tmp
            h1=  1.0+ 76.18009173 / (xx + 1)
            h2= -86.50532033 / (xx + 2)
            h3=  24.01409822 / (xx + 3)
            h4= -1.231739516 / (xx + 4)
            h5=  0.120858003e-2 / (xx + 5)
            h6= -0.536382e-5 / (xx + 6)
            ser= h1 + h2 + h3 + h4 + h5 + h6
            return tmp + np.log(2.50662827465 * ser)
        
        if self.s <= 0 or self.s >= 30 or not self.t.shape[0] > 1:
            hdir = np.zeros(np.shape(self.t))
            hdir[np.argmin(np.abs(self.t-self.t_mean))] = 1.
            self.dth = 1.
        else:
            help1 = np.exp((2 * self.s - 1) * np.log(2) +
                           2. * LnGamma(self.s + 1) -
                           LnGamma(2. * self.s + 1))
            help1  =help1 / np.pi
            help2 = np.abs(np.cos((self.t - self.t_mean) / 2)) ** (2 * self.s)
            hdir = help1 * help2

        if len(hdir) > 1:
            hdir /= np.sum(.5 * (hdir[1:] + hdir[:-1]) * self.dth)
        
        return np.outer(S, hdir)

        
if __name__ == "__main__":
    Nf = 10
    fstr = 0
    fstp = 0.5
    f = np.linspace(fstr,fstp,Nf)
    df = f[1]-f[0]
    fp = 1/9.
    Hs = 1.
    Wave_sp = wave_spec(f,fp,Hs)
    
    Jirr = 500.*Hs**2*(1./fp)
    Jreg = 1025*9.81**2/(32*np.pi)*Hs**2*(1/fp)
    
    #print('Energy flux for irregular waves')
    #print(Jirr)
    #print('*'*50)
    #print('Energy flux for regular waves')
    #print(Jreg)
    #print('*'*50)
    

    Wave_sp.spec_type = "Regular"
    Wave_sp.add_spectrum()
    ind = f>0
    J = 1025.*9.81**2/(4.*np.pi)*sum(Wave_sp.specs[-1][2][ind,0]/f[ind])*df
    #print(Wave_sp.specs_type_list[-1])
    #print(J)
    #print('*'*50)
    
    Wave_sp.spec_type = "Bretschneider_Mitsuyasu"
    Wave_sp.add_spectrum()
    J = 1025.*9.81**2/(4.*np.pi)*sum(Wave_sp.specs[-1][2][ind,0]/f[ind])*df
    #print(Wave_sp.specs_type_list[-1])
    #print(J)
    #print('*'*50)
    
    Wave_sp.spec_type = "Modified_Bretschneider_Mitsuyasu"
    Wave_sp.add_spectrum()
    J = 1025.*9.81**2/(4.*np.pi)*sum(Wave_sp.specs[-1][2][ind,0]/f[ind])*df
    #print(Wave_sp.specs_type_list[-1])
    #print(J)
    #print('*'*50)

    Wave_sp.spec_type = "Pierson_Moskowitz"
    Wave_sp.add_spectrum()
    J = 1025.*9.81**2/(4.*np.pi)*sum(Wave_sp.specs[-1][2][ind,0]/f[ind])*df
    #print(Wave_sp.specs_type_list[-1])
    #print(J)
    #print('*'*50)

    Wave_sp.spec_type = "Jonswap"
    Wave_sp.add_spectrum()
    J = 1025.*9.81**2/(4.*np.pi)*sum(Wave_sp.specs[-1][2][ind,0]/f[ind])*df
    #print(Wave_sp.specs_type_list[-1])
    #print(J)
    #print('*'*50)

    Wave_sp.spec_type = "pscTMA"
    Wave_sp.add_spectrum()

    Wave_sp.spec_type = "pscFRF"
    Wave_sp.add_spectrum()

    Wave_sp.spec_type = "pscSwell"
    Wave_sp.add_spectrum()
    
    Wave_sp.spec_type = "Regular"
    Wave_sp.add_spectrum()

    Wave_sp.show()
    
    Wave_sp.rm_spectrum(index='all')

    nt = 110
    Wave_sp.t = np.linspace(-np.pi,np.pi,nt)
    Wave_sp.spec_type = "Jonswap"
    Wave_sp.add_spectrum()
    
    s_vec = np.linspace(1,16,5)
    for s in s_vec:
        Wave_sp.s = s
        Wave_sp.add_spectrum()
        Wave_sp.show(index=[-1])
    

	
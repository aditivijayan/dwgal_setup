"""
Module for providing a galaxy+NFW+external gravity potential to the cooling_flow module
Vc modified with Agora vc
"""

import numpy as np
import scipy
from numpy import log as ln, log10 as log, e, pi, arange, zeros
from astropy import units as un, constants as cons
from astropy.cosmology import Planck15 as cosmo
import cooling_flow as CF
from scipy.interpolate import interp1d
from pathlib import Path

class NFW(CF.Potential):	
    mu=0.6
    X=0.75
    def __init__(self,Mvir,z,cvir,_fdr = 100.):
        self._fdr = _fdr
        self.Mvir = Mvir
        self.z = z
        self.cvir = cvir
        self.dr = self.r_scale()/self._fdr
        rs = arange(self.dr.value,self.rvir().value,self.dr.value) * un.kpc
        self.rho_scale = (self.Mvir / (4*pi * rs**2 * self.dr * 
                                               self.rho2rho_scale(rs) ).sum() ).to('g/cm**3') 

        vcirc = "vcirc_agora.txt"
        filename = Path(__file__).resolve().parent / vcirc
        data = np.loadtxt(filename)
        rad = data[:, 0]
        vcirc = data[:, 1]
        vcirc_interp = interp1d(rad,vcirc)
        Rmax = 5000. 
        r_grid = np.logspace(np.log10(0.001), np.log10(Rmax), 1000)
        vc_agora = np.zeros(r_grid.shape[0])
        i = 0
        for r in r_grid:
            if(r>np.amax(rad)):
                vc_agora[i] = vcirc[-1]
                
            elif(r<np.amin(rad)):
                vc_agora[i] = vcirc[0]
                
            else:
                vc_agora[i] = vcirc_interp(r)
            i+=1
        self.vc_interp = interp1d(r_grid,vc_agora)
            
    def Delta_c(self): #Bryan & Norman 98
        x = float(cosmo.Om(self.z) - 1)
        return float(18*pi**2 + 82*x - 39*x**2)
    def rvir(self):
        return ((self.Mvir / (4/3.*pi*self.Delta_c()*cosmo.critical_density(self.z)))**(1/3.)).to('kpc')
    def r_ta(self,use200m=False): 
        if not use200m:
            return 2*self.rvir()
        else:
            return 2*self.r200m()
    def r_scale(self):
        return self.rvir() / self.cvir
    def rho2rho_scale(self,r): 
        return 4. / ( (r/self.r_scale()) * (1+r/self.r_scale())**2 ) 
    def rho(self,r):
        return self.rho_scale * self.rho2rho_scale(r)
    def enclosedMass(self,r):
        return (16*pi*self.rho_scale * self.r_scale()**3 * 
                        (ln(1+r/self.r_scale()) - (self.r_scale()/r + 1.)**-1.)).to('Msun')
    def v_vir(self):
        return ((cons.G*self.Mvir / self.rvir())**0.5).to('km/s')
    def v_ff(self,r,rdrop=None):
        if rdrop==None: rdrop = 2*self.r200m()
        return ((2*(self.Phi(rdrop) - self.Phi(r)))**0.5).to('km/s')
    def vc(self,r):
        Ms = self.enclosedMass(r)
        vc = self.vc_interp(r) * un.km / un.s
        # return ((cons.G*Ms / r)**0.5).to('km/s')
        return vc
    def dlnvc_dlnR(self, r):
        return (ln(1+r/self.r_scale())*(self.r_scale()/r + 1.)**2 - (self.r_scale()/r + 1.))**-1.
    def mean_enclosed_rho2rhocrit(self,r):
        Ms = self.enclosedMass(r)
        return Ms / (4/3.*pi*r**3) / cosmo.critical_density(self.z)
    def r200(self,delta=200.):
        rs = arange(self.dr.value,2*self.rvir().value,self.dr.value)*un.kpc
        mean_rho2rhocrit = self.mean_enclosed_rho2rhocrit(rs)
        return rs[searchsorted(-mean_rho2rhocrit,-delta)]
    def r200m(self,delta=200.):
        rs = arange(self.dr.value,2*self.rvir().value,self.dr.value)*un.kpc
        mean_rho2rhocrit = self.mean_enclosed_rho2rhocrit(rs)
        return rs[searchsorted(-mean_rho2rhocrit,-delta*cosmo.Om(self.z))]		
    def M200(self,delta=200.):
        return self.enclosedMass(self.r200(delta))
    def M200m(self,delta=200.):
        return self.enclosedMass(self.r200m(delta))
    def Phi(self,r):
        return -(16*pi*cons.G*self.rho_scale*self.r_scale()**3 / r * ln(1+r/self.r_scale())).to('km**2/s**2')
    def g(self,r):
        Ms = self.enclosedMass(r)
        return cons.G*Ms / r**2
    def t_ff(self,r):
        return 2**0.5 * r / self.vc(r)

def searchsortedclosest(arr, val):
    if arr[0]<arr[1]:
        ind = np.searchsorted(arr,val)
        ind = minarray(ind, len(arr)-1)
        return maxarray(ind - (val - arr[maxarray(ind-1,0)] < arr[ind] - val),0)        
    else:
        ind = np.searchsorted(-arr,-val)
        ind = minarray(ind, len(arr)-1)
        return maxarray(ind - (-val + arr[maxarray(ind-1,0)] < -arr[ind] + val),0)        
def maxarray(arr, v):
    return arr + (arr<v)*(v-arr)
def minarray(arr, v):
    return arr + (arr>v)*(v-arr)
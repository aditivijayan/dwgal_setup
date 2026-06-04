"""
Module for providing the Wiersma et al. (2009) cooling functions to the cooling_flow module
"""
dataDir = 'cooling/'
import glob,h5py
import scipy, numpy as np
from scipy import integrate, interpolate
from numpy import log as ln, log10 as log, e, pi, arange, zeros
from astropy import units as un, constants as cons
import cooling_flow as CF
from pathlib import Path


class GrackleCooling(CF.Cooling):
    """
    creates Grackle cooling function for given metallicity and redshift
    """
    def __init__(self,Z2Zsun,z, shield=True):

        #Constants
        cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
        X = cloudy_H_mass_fraction
        Z = 0.02  # metal fraction by mass
        Y = 1. - X - Z
        mean_metals_A = 16.  # mean atomic weight of metals
        sigma_T = 6.6524e-25  
        electron_mass_cgs = 9.1093897e-28  # electron mass (g)
        mu = 0.6
        m_H = 1.672623e-24  # g
        T_cmb = 2.725  # * (1 + z); // K
        c_light_cgs_ = 2.99792458e10  # cgs
        radiation_constant_cgs_ = 7.5646e-15  # cgs
        boltzmann_constant_cgs_ = 1.380658e-16  # cgs
        
        if shield:
            table = 'CloudyData_UVB=HM2012_shielded.h5'
        else:
            table = 'CloudyData_UVB=HM2012.h5'
        filename = Path(__file__).resolve().parent / table
        with h5py.File(filename, "r") as f:
            metal_cooling = f["CoolingRates/Metals/Cooling"][:,0,:]
            metal_heating = f["CoolingRates/Metals/Heating"][:,0,:]
            primordial_cooling = f["CoolingRates/Primordial/Cooling"][:,0,:]
            primordial_heating = f["CoolingRates/Primordial/Heating"][:,0,:]

        net_metalprim_cooling_rate =  (Z2Zsun*(metal_cooling - metal_heating) + (primordial_cooling - primordial_heating))
        net_cooling_rate = np.transpose(net_metalprim_cooling_rate)
        nHbins = np.logspace(-10.,4, net_cooling_rate.shape[1])
        Tbins  = np.logspace(1.,9, net_cooling_rate.shape[0])
        Tgrid, nHgrid = np.meshgrid(Tbins, nHbins, indexing='ij')

        rhoHgrid = nHgrid* m_H
        rho = rhoHgrid / cloudy_H_mass_fraction


        #Photoelectric Cooling Rate
        n_e = (rho / m_H) * \
                         (1.0 - mu * (X + Y / 4. + Z / mean_metals_A)) / \
                         (mu - ( electron_mass_cgs / m_H))
        
        Tsqrt = np.sqrt(Tgrid)
        phi = 0.5  # phi_PAH from Wolfire et al. (2003)
        G_0 = 1.7  # ISRF from Wolfire et al. (2003)
        epsilon = \
        4.9e-2 / (1. + 4.0e-3 * (G_0 * Tsqrt / (n_e * phi))**0.73) + \
        3.7e-2 * (Tgrid / 1.0e4)**(0.7) / \
            (1. + 2.0e-4 * (G_0 * Tsqrt / (n_e * phi)))
        Gamma_pe = 1.3e-24 * nHgrid * epsilon * G_0
        cooling_rate_PE = -1. * (Z2Zsun * Gamma_pe)  / nHgrid ** 2

        #Compton Cooling
        E_cmb = radiation_constant_cgs_ * (T_cmb * T_cmb * T_cmb * T_cmb)
        Gamma_C = (8. * sigma_T * E_cmb) / (3. * cons.m_e.value * c_light_cgs_)
        C_n     = Gamma_C * boltzmann_constant_cgs_ / (5. / 3. - 1.0)
        compton_CMB =  (C_n * (Tgrid - T_cmb) * n_e / nHgrid ** 2)

        net_lambda = (net_cooling_rate + cooling_rate_PE + compton_CMB)

                
        self.f_Cooling = interpolate.RegularGridInterpolator((log(Tbins), log(nHbins)),
                                                        net_lambda, 
                                                        bounds_error=False, fill_value=None)
        #### calculate gradients of cooling function
        X, Y = np.meshgrid(Tbins, nHbins, copy=False)
        # dlogT = np.diff(log(Tbins))[0] 
        # dlogn = np.diff(log(nHbins))[0] 
        # vals = log(self.LAMBDA(X*un.K,Y*un.cm**-3).value)
        # dLambda_drho, dLambda_dT = np.gradient(vals, nHbins, Tbins)
        # dlnLambda_dlnrhoArr, dlnLambda_dlnTArr = np.gradient(vals,dlogn, dlogT)  
        logTbins  = np.log(Tbins)
        lognHbins = np.log(nHbins)
        dLambda_dlnT, dLambda_dlnnH = np.gradient(net_lambda, logTbins, lognHbins)
        dlnLambda_dlnrhoArr = (1./net_lambda) * dLambda_dlnnH
        dlnLambda_dlnTArr   = (1./net_lambda) * dLambda_dlnT
        
        self.dlnLambda_dlnT_interpolation = interpolate.RegularGridInterpolator((log(Tbins), log(nHbins)),dlnLambda_dlnTArr, bounds_error=False, fill_value=None)
        self.dlnLambda_dlnrho_interpolation = interpolate.RegularGridInterpolator((log(Tbins), log(nHbins)),dlnLambda_dlnrhoArr, bounds_error=False, fill_value=None)                
    def LAMBDA(self, T, nH):
        """cooling function"""
        return self.f_Cooling((log(T.to('K').value), log(nH.to('cm**-3').value))) * un.erg*un.cm**3/un.s
    def tcool(self,T,nH):
        """cooling time"""
        return 3.5 * cons.k_B * T / (nH * self.LAMBDA(T, nH))
    def f_dlnLambda_dlnT(self, T, nH):         
        """logarithmic derivative of cooling function with respect to T"""
        return self.dlnLambda_dlnT_interpolation((log(T.to('K').value), log(nH.to('cm**-3').value)))
    def f_dlnLambda_dlnrho(self, T, nH):
        """logarithmic derivative of cooling function with respect to rho"""
        return self.dlnLambda_dlnrho_interpolation((log(T.to('K').value), log(nH.to('cm**-3').value)))
    



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


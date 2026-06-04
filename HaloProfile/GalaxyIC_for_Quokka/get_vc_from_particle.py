import sys
from pathlib import Path
sys.path.append(str(Path.cwd().parent))
import os
import pandas as pd

from astropy import units as u, constants as const
import numpy as np
import HaloPotential as Halo
from scipy.interpolate import interp1d

Mvir=1.775e11*const.M_sun #LMC Mass from Shah+
potential = Halo.NFW(Mvir=Mvir, z=0, cvir=10)
Rvir = potential.rvir().value   

G = const.G.cgs.value 
kpc = (1 * u.kpc).to(u.cm).value
kmps = (1 * u.km / u.s).to(u.cm / u.s).value 


def acc_due_to_grav(x, x_part, y_part, z_part, mass):
    denom = ((x-x_part)**2 + y_part**2 + z_part**2)**1.5
    numer = G * mass * (x_part - x)
    acc_g = np.sum(numer/denom)
    return acc_g
    

def make_vc_interp(infile):
    x_part,y_part,z_part,mass = np.loadtxt(infile, dtype=float, comments='#', delimiter=None, skiprows=1, usecols=(0, 1, 2,3), unpack=True)
    x = np.logspace(-1, 4, 100) * kpc
    g = np.zeros(x.shape[0])
    i=0
    for x1 in x:
        g[i] = acc_due_to_grav(x1, x_part,y_part,z_part,mass)
        i+=1
    vc2 = np.abs(g) * x
    vc = vc2**0.5 / kmps
    window_size = 5
    window = np.ones(window_size) / window_size
    vc_smoothed = np.convolve(vc, window, mode="same")
    log_vc_interp = interp1d(np.log10(x), np.log10(vc_smoothed), kind="cubic", bounds_error=False, fill_value="extrapolate")
    return log_vc_interp

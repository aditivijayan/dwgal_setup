
import sys
from pathlib import Path
sys.path.append(str(Path.cwd().parent))


from astropy import units as un, constants as cons
import numpy as np
import matplotlib.pyplot as plt
from constants import *


def set_plot_style():
    plt.rcParams.update({
        "font.size": 22,
        "axes.linewidth": 2,

        "xtick.major.size": 10,
        "xtick.minor.size": 5,
        "xtick.major.width": 2,
        "xtick.minor.width": 1,
        "xtick.direction": "in",

        "ytick.major.size": 10,
        "ytick.minor.size": 5,
        "ytick.major.width": 2,
        "ytick.minor.width": 1,
        "ytick.direction": "in",
    })

def make_plots(prof, anl_profile):
    
    Mdot = prof["Mdot"]
    Rvir = prof["Rvir"]
    Mvir = prof["Mvir"]
    
    fig, ax = plt.subplots(4, 1, gridspec_kw = {'wspace':0.1, 'hspace':0.0},figsize=(6, 18))

    ylabel = [r'$\rho/m_p$', r'$T$ [K]', r'$M$', r'$t_{\rm cool}$ [Gyr]']
    
    ax[0].plot(prof["r"], prof["rho"]*X/mp)
    ax[0].plot(prof["r"], anl_profile["nH"])

    ax[1].plot(prof["r"], prof["T"], label='Solution')
    ax[1].plot(prof["r"], anl_profile["T"], label='Analytical')
    
    ax[2].plot(prof["r"], prof["Mach"])
    ax[2].plot(prof["r"], anl_profile["Mach"])
    
    ax[3].plot(prof["r"], prof["tcool"])
    ax[3].plot(prof["r"], anl_profile["tcool"])
    
    
    ax[-1].set_xlabel('r [kpc]')
    for i in range(4):
        ax[i].set_ylabel(ylabel[i])
    plt.setp(ax, 'xscale', ('log'))
    plt.setp(ax, 'yscale', ('log'))
    
    ax[0].set_ylim(1.e-7, 2.e-3)
    ax[1].set_ylim(2.e4, 8.e6)
    ax[2].set_ylim(1.e-2,1.e1)
    ax[3].set_ylim(1.e-3,90.)
    
    ax[0].tick_params(axis='x', which='both', labelbottom=False, top=True, bottom=True)
    ax[1].tick_params(axis='x', which='both', labelbottom=False, top=True, bottom=True)
    ax[-2].tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=False)
    ax[-1].tick_params(axis='x', which='both', top=True, bottom=True, labelbottom=True)
    ax[0].set_title(r'$\dot{M}$=%.1e'%(Mdot*yr_to_sec/Msun) + r' $M_{\odot} \ \rm{yr}^{-1}$')
    ax[0].text(0.4, 0.85, r'$M_{\rm halo}$=%.3e'%(Mvir/Msun), transform=ax[0].transAxes)
    ax[0].text(0.4, 0.75, r'$r_{\rm min}$=%.1f'%(np.amin(prof["r"])), transform=ax[0].transAxes)
    ax[0].text(0.4, 0.65, r'$R_{\rm vir}$=%.1f'%(Rvir) + ' kpc', transform=ax[0].transAxes)
    ax[1].legend()
    return fig

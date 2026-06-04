import sys
from pathlib import Path
sys.path.append(str(Path.cwd().parent))
import argparse

from astropy import units as un, constants as cons
import numpy as np
import cooling_flow as CF
import HaloPotential as Halo
import matplotlib.pyplot as plt
import GrackleCooling as Cool
import AgoraPotential as Agora
from constants import *
from make_plots import make_plots, set_plot_style
from scipy.interpolate import RegularGridInterpolator
from get_vc_from_particle import make_vc_interp


def setup_halo_and_cooling(Mvir, cvir, Z2Zsun, z, shield):
    potential = Halo.NFW(Mvir=Mvir, z=z, cvir=cvir)

    Z2Zsun = Z2Zsun
    z = z
    cooling = Cool.GrackleCooling(Z2Zsun, z, shield=shield)

    return potential, cooling


def solve_transonic_flow(potential, cooling, R_sonic, R_min=None, R_max=None, max_step=0.1):
    solution = []
    if R_max is None:
        R_max = potential.rvir()

    solution.append( CF.shoot_from_sonic_point(
        potential,
        cooling,
        R_sonic,
        R_max,
        R_min,
        max_step=max_step,
        calcInwardSolution=False,
        pr=True
    ))
    if(len(solution) == 0):
        raise RuntimeError("No solution found. Check parameters.")
    else:
        res = solution[0]
        return res


def compute_profiles(res, potential):
    r = res.Rs().value
    rho = res.rhos().value
    cs = res.cs().value
    Mdot = res.Mdot.value * (Msun / yr_to_sec)
    v = Mdot / (4*np.pi*(r*kpc)**2 * rho)
    Mach = v / cs / kmps
   
    T = cs**2 * mu * mp * kmps**2 / (gamma * kb)
    tcool = res.t_cools().value
    lambda_cool = np.abs(res.Lambdas().value)

    path_to_particle_file = '/Users/aditivijayan/quokka/INCITE/dwGal_setup/GalIC'
    file_name = 'LMC_1000sol.txt'
    infile = os.path.join(path_to_particle_file, file_name)
    log_vc_interp = make_vc_interp(infile)
    vc = 10.**log_vc_interp(np.log10(r*kpc))
    return {
        "r": r,
        "rho": rho,
        "T": T,
        "v": v,
        "vc": vc,
        "Mach": Mach,
        "tcool": tcool,
        "lambda": lambda_cool,
        "Mdot": Mdot,
        "Rvir": potential.rvir().value,
        "Mvir": potential.Mvir.to("g").value
    }


def analytic_profiles(r, Mhalo, Mdot, lambda_cool):
    Mhalo = Mhalo.to("g").value
    Mach_anl = 0.11 * (Mhalo / (1e12 * Msun))**(-0.72) \
                     * (Mdot / (Msun / yr_to_sec))**0.5 \
                     * (lambda_cool / 1e-22)**0.5 \
                     * (r / 100)**(-0.3)

    nH_anl = 1.6e-5 * (Mhalo / (1e12 * Msun))**0.36 \
                     * (Mdot / (Msun / yr_to_sec))**0.5 \
                     * (lambda_cool / 1e-22)**(-0.5) \
                     * (r / 100)**(-1.6)

    tcool_anl = 7.2 * (Mhalo / (1.e12 * Msun))**0.36 \
                    * (Mdot / (Msun / yr_to_sec))**(-0.5) \
                    * (lambda_cool / 1.e-22)**(-0.5) \
                    * (r / 100)**1.4

    T_anl = (
        X * mu * tcool_anl * (gamma - 1)
        * nH_anl * lambda_cool * 1.e9 * yr_to_sec / kb
    )

    analy_profile = {
        "Mach": Mach_anl,
        "nH": nH_anl,
        "T": T_anl,
        "tcool": tcool_anl,
    }

    return analy_profile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Solve transonic cooling-flow solution"
    )

    parser.add_argument(
        "--Mvir",
        type=float,
        default=1.775e11,
        help="Halo virial mass in Msun"
    )

    parser.add_argument(
        "--c",
        type=float,
        default=10.0,
        help="NFW concentration parameter"
    )

    parser.add_argument(
        "--Z2Zsun",
        type=float,
        default=0.5,
        help="Metallicity relative to solar"
    )

    parser.add_argument(
        "--Rsonic",
        type=float,
        default=2.9,
        help="Sonic radius in kpc"
    )

    parser.add_argument(
        "--shield",
        type=bool,
        default=True,
        help="If True use shielding in cooling function"
    )

    parser.add_argument(
        "--z",
        type=float,
        default=0.0,
        help="Redshift"
    )

    return parser.parse_args()

def write_ic_file(filename, prof, potential, Z2Zsun, rsonic):
    data = np.column_stack([
        prof["r"],
        prof["vc"],
        prof["rho"],
        prof["v"],
        prof["T"],
    ])

    header = ("r[kpc] vc[km/s] rho v T" 
        f"Mdot[Msun/yr]={prof['Mdot']*yr_to_sec/Msun:.6e} "
        f"rsonic[kpc]={rsonic:.6e} "
        f"Z2Zsun={Z2Zsun:.2e}"
    )
    
    np.savetxt(filename, data, header=header, fmt="%.6e")

def main():
    
    args = parse_args()

    # User inputs → physical quantities
    Mvir = args.Mvir * cons.M_sun
    cvir = args.c
    Z2Zsun = args.Z2Zsun
    Rsonic = args.Rsonic * un.kpc
    z = args.z
    shield = args.shield
    print(
    f"Generating ICs for Mvir={Mvir/cons.M_sun:.3e}, "
    f"c={cvir}, Z2Zsun={Z2Zsun}, Rsonic={Rsonic}, z={z}"
)


    # Setup potential
    potential, cooling = setup_halo_and_cooling(Mvir=Mvir, cvir=cvir, Z2Zsun=Z2Zsun, z=z, shield=shield)

   #Solve for cooling flows
    res = solve_transonic_flow(potential, cooling, R_sonic= Rsonic)
    profiles = compute_profiles(res, potential)

    anl_profiles = analytic_profiles(
        profiles["r"],
        Mvir,
        profiles["Mdot"],
        profiles["lambda"]
    )


    #Plot solutions and save data
    set_plot_style()

    fig = make_plots(profiles, anl_profiles)

    image_name = "solution_%.2e" % (profiles["Mdot"]*yr_to_sec/Msun) + ".jpeg" 
    data_file = "lmc_halo_IC.dat"
    print("Saving figure... in", image_name)
    fig.savefig(image_name, dpi=300, bbox_inches='tight')
    print("Writing IC file... in", data_file)
    write_ic_file(data_file, profiles, potential, Z2Zsun, Rsonic.value)


if __name__ == "__main__":
    main()
    

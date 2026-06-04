import numpy as np
import matplotlib.pyplot as plt
import os,re
import h5py, argparse
from readparam import read_galic_params
from pathlib import Path

Msun = 1.989e+33
kmps = 1.e5
kpc  = 3.018e21


parser = argparse.ArgumentParser(description="Process snapshot folder")

parser.add_argument(
    "--boxsize",
    type=float,
    required=True,
    help="Box size (in kpc)"
)

parser.add_argument(
    "--param",
    type=str,
    required=True,
    help="Name of the parameter file"
)

args = parser.parse_args()

param_file   = Path(args.param)
boxsize = args.boxsize

output_dir_from_param, mass_unit = read_galic_params(param_file)
input_folder = (param_file.parent / output_dir_from_param).resolve()


snap_re = re.compile(r"snap_(\d+)\.hdf5$")

max_snap = -1
latest_file = None

for fname in os.listdir(input_folder):
    m = snap_re.match(fname)
    if m:
        snap_num = int(m.group(1))
        if snap_num > max_snap:
            max_snap = snap_num
            latest_file = fname

if latest_file is None:
    raise RuntimeError("No snap_XX.hdf5 files found")

infile = os.path.join(input_folder, latest_file)

print("Analysing snapshot:", infile)

hf = h5py.File(infile, 'r')

#Read Halo particles
coord_halo = hf["PartType1"]["Coordinates"][:]
x_halo = coord_halo[:,0]
y_halo = coord_halo[:,1]
z_halo = coord_halo[:,2]
plt.scatter(x_halo, y_halo)
plt.xlim(-boxsize, boxsize)
plt.ylim(-boxsize, boxsize)

vel_halo = hf["PartType1"]["Velocities"][:]
vx_halo = vel_halo[:,0]
vy_halo = vel_halo[:,1]
vz_halo = vel_halo[:,2]

mass_halo = hf["PartType1"]["Masses"][:]
mass_halo_part = mass_halo.astype(np.float64) * mass_unit
print("Mass of halo particles:", np.amax(mass_halo_part)/Msun, ", ", np.amin(mass_halo_part)/Msun)

#Read Disc particles
coord_disc = hf["PartType2"]["Coordinates"][:]
x_disc = coord_disc[:,0]
y_disc = coord_disc[:,1]
z_disc = coord_disc[:,2]
plt.scatter(x_disc, y_disc)
plt.xlim(-boxsize, boxsize)
plt.ylim(-boxsize, boxsize)

vel_disc = hf["PartType2"]["Velocities"][:]
vx_disc = vel_disc[:,0]
vy_disc = vel_disc[:,1]
vz_disc = vel_disc[:,2]

mass_disc = hf["PartType2"]["Masses"][:]
mass_disc_part = mass_disc.astype(np.float64) * mass_unit
print("Mass of disc particles:", np.amax(mass_disc_part)/Msun, ", ", np.amin(mass_disc_part)/Msun)

x_all = np.concatenate([x_halo, x_disc])
y_all = np.concatenate([y_halo, y_disc])
z_all = np.concatenate([z_halo, z_disc])
vx_all = np.concatenate([vx_halo, vx_disc])
vy_all = np.concatenate([vy_halo, vy_disc])
vz_all = np.concatenate([vz_halo, vz_disc])
mass_all = np.concatenate([mass_halo_part, mass_disc_part])


# Create a boolean mask: True for particles inside the box
mask = (x_all > -boxsize) & (x_all < boxsize) & \
       (y_all > -boxsize) & (y_all < boxsize) & \
       (z_all > -boxsize) & (z_all < boxsize)

# Apply the mask
x_box = x_all[mask]
y_box = y_all[mask]
z_box = z_all[mask]

vx_box = vx_all[mask]
vy_box = vy_all[mask]
vz_box = vz_all[mask]

mass_box = mass_all[mask]
print(f"Selected {len(x_box)} particles out of {len(x_all)} total particles")
filename = f"{output_dir_from_param}.txt"
print("Saving to file:", filename)


header = str(x_box.shape[0])
data = np.column_stack([x_box*kpc, y_box*kpc, z_box*kpc, mass_box, vx_box*kmps, vy_box*kmps, vz_box*kmps ])
np.savetxt(filename, data, header=header, fmt="%.6e", comments='')

# dwGal_setup

This repository contains tools and configurations to set up initial conditions (ICs) for Dwarf Galaxy simulations within the **Quokka** hydrodynamics code.

The repository is organized into two main components:

---

## 📁 Repository Structure

### 1. HaloProfile
This component contains the cooling flow solutions based on the **Stern et al. (2019)** models. 
* **Primary Use:** We utilize these solutions **strictly for establishing the gravitational potential** of the host halo.
* **CGM Physics:** For specific dwarf galaxy runs (such as the Large Magellanic Cloud `LMC` or `M82`), we assume a **hydrostatic, isothermal Circumgalactic Medium (CGM)** fixed at a temperature of $3 \times 10^6\text{ K}$ ($3\text{ MK}$).

### 2. GalIC_for_Quokka
This directory contains the pipeline and processing scripts needed to convert raw particle data into Quokka-compatible initial condition grids.
* **Functionality:** It processes standard **GalIC** particle files and formats them into the specific grid layouts or HDF5 files required by Quokka’s setup routines.
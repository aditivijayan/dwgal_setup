# This folder generates the [Stern+19](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2549S/abstract) cooling flow solutions. Cooling flow solutions can be generated given a galaxy halo mass.
# I have added the following functionalities on top of the original code in order to generate halo profiles in a Quokka-compatible format:

    1. AgoraPotential.py : This file contains a class that uses vcirc_agora.txt to create an Agora-like potential.
    2. solve_for_profiles.py : This file contains the main driver code to generate cooling flow solutions.
        -- Halo mass, concentration, metallicity, sonic radius and redshift are passed by the user. The default values correspond to LMC-like halo.
        -- The code then uses GrackleCooling.py to create a cooling function with/without shielding based on user input.
        -- It calls the cooling flow solver to generate the cooling flow solution.
        -- Finally, it saves the profiles in a Quokka-compatible format


# Things to note- 
     -- Make sure the halo mass are passed correctly.
     -- Rmax is taken to be Rvir.
     -- If no solutions are found for given sonic radius, try changing the sonic radius. I have found that it is 2.9 kpc for LMC and 0.3 kpc for MW halo.
     -- The solution is not valid for all r. It holds only for Rsonic < r < Rvir.       
     -- By default no solutions are found for r< Rsonic. But you can change this setting by calcInwardSolution=True in solve_transonic_flow function.
     -- The analytical solutions are taken from Stern+19. 
     -- I personally have not tested solve_profile.py. But I have tested the module using hydro-static-isothermal-CGM.ipynb.

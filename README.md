# A parameter representing missing charge should be considered when calibrating action potential models

This GitHub repo was design to facilitate the reproduction of the results published in the paper linked to this GitHub.

Note that to be able to run the Python code, one needs to have installed the Python packages:
  - Myokit v1.30.2: https://github.com/MichaelClerx/myokit 
  Myokit: A simple interface to cardiac cellular electrophysiology Michael Clerx, Pieter Collins, Enno de Lange, Paul   G.A. Volders 2016 Progress in Biophysics and Molecular Biology Volume 120, issues 1-3, pages 100-114, 
  doi:  10.1016/j.pbiomolbio.2015.12.008
  
  - SABS_project: https://github.com/rcw5890/SABS_project
  This package is an interface to Myokit that was used for the present project. Note that Myokit is the engine that sets up the simulation, and solves it using SUNDIALS solver.

  - Matplotlib v3.5.1. 
  Use the classic ```pip install --upgrade matplotlib``` command to install the latest version.

All the scripts that were used to generate the figures in the present article can be found in the ```Scripts``` folder. To reproduce the figures, 1) download the repository, 2) in the scripts, change the working directory to the folder where the repository is located, 3) run the code with ```python 'name_of_the_script.py'```!

## Note concerning the reproduction of the fitting of the ORd-CiPA model

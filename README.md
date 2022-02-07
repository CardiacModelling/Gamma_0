# A parameter representing missing charge should be considered when calibrating action potential models

This GitHub repo was design to facilitate the reproduction of the results published in the paper linked to this GitHub.

## List of Python packages necessary for the reproduction of this work
To be able to run the Python scripts, one needs to have installed:
  - Python v3.7.4 or later version 
  - Myokit v1.30.2: https://github.com/MichaelClerx/myokit 
  Myokit: A simple interface to cardiac cellular electrophysiology Michael Clerx, Pieter Collins, Enno de Lange, Paul   G.A. Volders 2016 Progress in Biophysics and Molecular Biology Volume 120, issues 1-3, pages 100-114, 
  doi:  10.1016/j.pbiomolbio.2015.12.008
 
  - SABS_project: https://github.com/rcw5890/SABS_project
  This package is an interface to Myokit that was used for the present project. Note that Myokit is the engine that sets up the simulation, and solves it using SUNDIALS solver.
  - Matplotlib v3.5.1. 
  Use the classic ```pip install --upgrade matplotlib``` command to install the latest version.
  
All the scripts that were used to generate the figures in the present article can be found in the ```Scripts``` folder. To reproduce the figures, 1) download the repository, 2) in the scripts, change the working directory to the folder where the repository is located, 3) run the code !

## Note for new Python users
For users who are new to Python, here are more details on how to install the Python packages and run the scripts.

1) Install Anaconda: https://www.anaconda.com/products/individual
2) In the installation process, install Spyder.
3) Download the needed packages.
4) Open the Anaconda prompt to install the packages.
  - Update pip: use the command ```pip install --upgrade pip --user```. The ```--user``` argument might not be needed if you have the admin rights for your session.
  - Navigate using the ```cd your_directory``` command to get to the location of the ```setup.py``` script of the downloaded library.
  - Install the package with the command ```pip install . -e --user```. Note that for the installation of Myokit, you need to have a C++ compiler installed. More details concerning the installation of Myokit can be found at http://myokit.org/install.
5) In the Spyder interpreter, load the script of interest.
6) In the first compartment, change the line starting with ```os.chdir(``` to change the folder name.
7) Run the script!

## Note concerning the reproduction of the fitting of the ORd-CiPA model
The scripts to fit the ORd-CiPA model were executed on a HPC owned by Roche, with parallelisation on 7 cores. The execution took between one and two hours. The CSV files outputs are included in this repository.
To reproduce the results on a Windows machine, it can be necessary to turn off the parallel evaluation of proposed parameters. Set ```parallel = False``` at the beginnning of the script in that case.

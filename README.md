# ISW Lab Toolkit

This is a toolkit/collection of scripts, functions, bits of files that help to make the lab work go smoothly in the Geophysical Fluid Dynamics lab in Newcastle University.

## Introduction
The lab data primarily comes in 2 forms, the recorded videos from cameras, processed initially in DigiFlow, and the density data from microconductivity probes (which run on an Arduino system). The scripts contain the main functions to read this data into MATLAB and analyse, as well as in places scripts that interface into other software (primarily here DigiFlow). A couple of the MATLAB scripts help in the basic day-to-day running of the lab, e.g. `copycamfiles.m`, `Camera_parallelisation.m` 

## Installation 
### Full Installation
First off, 
	```git clone --recurse-submodules https://github.com/HartharnSam/ISWLabToolkit```
Then you'll want to add this whole ordeal to your MATLAB path, which can be done:
- temporarily with  ```addpath(genpath('path/to/the/toolkit'));```
- permanently by adding the addpath... code to `startup.m`. startup.m should be found in the file you find by typing userpath into the matlab command window - if it is not, create one.

## Files:
- /Cameras *This directory should contain files needed for camera operation in DigiFlow. More detailed README within*
- /DigiFlow_read *This directory should contain files needed for reading & processing DigiFlow data in MATLAB. More detailed README within*
- /DJLES *DJL Solver for MATLAB. Submodule only - so links directly to the github repository.* 
- /Probes *More detailed README within*
	- /Probe Drivers *This directory should contain files needed for operating the probes using the arduino software. 
	- Main directory contains files needed for processing probe data in MATLAB
- `calc_DJL.m` *Calculates a DJL wave to match a given lab setup (e.g. actual densities, layer depths & wave amplitude)*
- `calc_kdv.m` *Same as calc_DJL but for a kdv wave*
- `Camera_parallelisation.m` *helps in the task of parallelising the cameras to the tank*
- `copycamfiles.m` *copies the day's data from local D: drives to a central OneDrive - always needs modifications to new systems*
- `LabWaveProperties.m` *New method of picking points on the pycnocline using MATLAB GUI, and saves/outputs amplitudes and wave velocity. Crucially, saves the picked points for replicability*
- `setup.sh` *Sets up folder structures each day, makes sure copies of relevant scripts are where they'll need to be later on*

## Wave Solutions

Two scripts are included here (plus the entire DJLES package) which provide analytical solutions to ISWs. This is the Weakly Nonlinear Korteweg de Vries theory, and the fully nonlinear Dubreil-Jacotin-Long (DJL) equation.

### KdV
We get:

$$ \lambda(-D^2+k^2 I)\phi = N^2(z)\phi $$

$$ B_t = -cB_x + \epsilon 2r_{10}cBB_x + \mu r_{01}B_{xxx} $$

Where $\phi(z)$ is the vertical structure function, $\epsilon$ is a small dimensionless parameter measuring wave amplitude, $\mu = (H/L)^2$ is a measure of the aspect ratio.

Further discussion on the KdV equation can be found in [Lamb \& Yan 1996](https://doi.org/10.1175/1520-0485(1996)026%3C2712:TEOIWU%3E2.0.CO;2). 


### DJL
Fully nonlinear equation from the steady, incompressible Euler equations. In a constant background current:

$$\nabla^2 \eta + \frac{N^2(z-\eta)}{U_0^2}\eta = 0$$

## Acknowledgements & Licence
This really is a collection of projects and scripts from various authors. The project overall was created by Sam Hartharn-Evans. The initial dfireadvel package, which became digiflow_read was created by J.K.Sveen@damtp.cam.ac.uk. The initial dissipation package within DigiFlow_read was created by Peter Sutherland. There are also contributions from Marek Stastna (University of Waterloo), Magda Carr (Newcastle University) and the DJLES package is the Dunphy et al. (2011) package. 
Where possible, individual acknowledgements have been made within scripts. 

The package is licenced under the MIT licence, see `LICENSE` for more information. 


## Further Reading
Things I found useful in this code writing:
- [MatPIV](https://www.mn.uio.no/math/english/people/aca/jks/matpiv/): A package which does the same PIV and PTV calculations in MATLAB. Helped me understand what's going on.
- [ISW_ParticleTrackModel](https://github.com/HartharnSam/ISW_ParticleTrackModel) : A Package for tracking particles at the surface as part of the ISW in Ice Covered Waters Project (ONE Planet PhD Studentship)

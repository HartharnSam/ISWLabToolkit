# ISW Lab Toolkit

This is a toolkit/collection of scripts, functions, bits of files that help to make the lab work go smoothly in the Geophysical Fluid Dynamics lab in Newcastle University. The bulk of most subpackages have been written by other authors and then re-purposed here, acknowledgements are found within the subpackage readmes and individual scripts where possible. 

## Introduction
The lab data primarily comes in 2 forms, the recorded videos from cameras, processed initially in DigiFlow, and the density data from microconductivity probes (which run on an Arduino system). The scripts contain the main functions to read this data into MATLAB and analyse, as well as in places scripts that interface into other software (primarily here DigiFlow). A couple of the MATLAB scripts help in the basic day-to-day running of the lab, e.g. `copycamfiles.m`, `Camera_parallelisation.m` 

There are some other files that are useful, e.g. parameter calculating spreadsheets and records sheets. TODO: Upload these. 

## Installation 
### Full Installation
`git clone --recurse-submodules https://github.com/HartharnSam/ISWLabToolkit`
Then you'll want to add this whole ordeal to your MATLAB path (temporarily with `addpath(genpath('path/to/the/toolkit'));`) or permanently by adding that code to `startup.m`

### Partial Installation 
TODO

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

## Acknowledgements & Licence
This really is a collection of projects and scripts from various authors, where possible acknowledgements have been made within scripts. 

The overall package is licenced as MIT. 


## Further Reading
Things I found useful in this code writing:
- [MatPIV](https://www.mn.uio.no/math/english/people/aca/jks/matpiv/): A package which does the same PIV and PTV calculations in MATLAB. Helped me understand what's going on.
- [ISW_ParticleTrackModel](https://github.com/HartharnSam/ISW_ParticleTrackModel) : A Package for tracking particles at the surface as part of the ISW in Ice Covered Waters Project (ONE Planet PhD Studentship)

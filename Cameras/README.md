# Cameras
This directory should contain files needed for camera operation in DigiFlow.

## Files
- [UniqVision_UP1830CL_8bit.r64](./UniqVision_UP1830CL_8bit.r64) & [UniqVision_UP1830CL_12bit.r64](./UniqVision_UP1830CL_8bit.r64): Camera configuration files for Uniq Vision UP1830 CL cameras in 8bit and 12bit mode. 8bit version should be used. A copy of this file should be made into `C:/BitFlow SDK 6.40/Config/R64/`
- [twocameras.dfc](./twocameras.dfc):  DigiFlow Code to simultaneously operate two cameras (on the same BitFlow board, different ports). 
- [dft2txt.dfc](./dft2txt.dfc): DigiFlow Code to export .dft format files (produced in Particle Tracking Velocimetry) to more reader-friendly .txt formats

## Dependencies
- A current installation guide for DigiFlow \& BitFlow can be found [here](https://github.com/HartharnSam/DigiFlow-Installation) , but note this is in a private repository, due to links to restricted links - please contact the owner of this repository to get access! 
- We use BitFlow SDK6.40. Newer versions may exist, but we know this version works

## twocameras
### Setup
Move the .dfc file to the parent directory for the current experiment directory. First time setup will require creation of a folder V:\Capture, with Full Control for Authenticated Users
### Usage
Check the board numbers are set up as intended. Open the dfcConsole (Edit>dfcConsole, or ctrl+E), and open the twocameras.dfc file. The code is inexplicably glitchy, and needs to be run section at a time, with sections indicated by #### (indicated by ##### markers). Select the section to be run, and Execute Section 

## dft2txt
This code runs through .dft files produced for each frame in Particle Tracking Velocimetry, it opens them and re-formats the data so it can be stored in two text files per frame - ptv_basic_####.txt and ptv_full_#####.txt. 

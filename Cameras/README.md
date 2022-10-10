This directory should contain files needed for camera operation in DigiFlow.

## Files
- UniqVision_UP1830CL_8bit.r64 & UniqVision_UP1830CL_12bit.r64 : Camera configuration files for Uniq Vision UP1830 CL cameras in 8bit and 12bit mode. 8bit version should be used. A copy of this file should be made into C:/
- twocameras.dfc :  DigiFlow Code to simultaneously operate two cameras (on the same BitFlow board, different ports). Note the code is glitchy, and needs to be run section at a time (indicated by ##### markers)

- Current copy of DigiFlow can be found at:
	https://www.dropbox.com/sh/uc5nmllmd11pdj1/AADNkwK5MmsYD2JW1gcFcXbTa?dl=0
- We use BitFlow SDK6.40. Newer versions may exist, but we know this version works


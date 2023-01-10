Scripts related to the operation of conductivity probes through the Arduino system. 
Contributions from Paul Watson (Newcastle University)

## Files:
- Probe Drivers: arduino board code for collecting dataManualProbeOperation: Manually operated probes, press button for each batch of readings at a given and manually recorded height
	- [`ProbeOperation.ino`](<./Probe Drivers/ProbeOperation/ProbeOperation.ino>): For normal probe operation with functioning potentiometer
	- [`TestPotentiometer.ino`](<./Probe Drivers/TestPotentiometer/TestPotentiometer.ino>): Prints outputs every second for general diagnostics
	- [`ManualProbeOperation.ino`](<./Probe Drivers/ManualProbeOperation/ManualProbeOperation.ino>): Returns outputs for each button press, with manually recorded depths
- [`probe_read.m`](./probe_read.m) : Matlab function for reading in probe data
- [`ProbeParse.m`](./ProbeParse.m) : Matlab script placed inside folder containing each probe's data which parses data to probe_read.m to output a _.mat_ file containing a structure relateing to the density profile; plot the profiles; and print the fitting parameters
- [`calc_N2.m`](./calc_N2.m) : Calculate the $N^2$ buoyancy frequency profile for a given density profile
	

	
	

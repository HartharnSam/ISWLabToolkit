Scripts related to the operation of probes through the Arduino system. 
Contributions from Paul Watson (Newcastle University)

## Files:
- Probe Drivers: arduino board code for collecting dataManualProbeOperation: Manually operated probes, press button for each batch of readings at a given and manually recorded height
	- [`ProbeOperation.ino`](./Probe Drivers/ProbeOperation/ProbeOperation.ino): For normal probe operation with functioning potentiometer
	- PotentiometerTest: Prints outputs every second for general diagnostics
- probe_read.m : Matlab function for reading in probe data
- ProbeParse.m : Matlab script placed inside folder containing each probe's data which parses data to probe_read.m to output a density profile structure and plot
	

	
	

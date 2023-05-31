# Probes
Scripts related to the operation of conductivity probes through the Arduino system. 
Contributions from Paul Watson (Newcastle University)

## Files:
### Probe Drivers
Each file for running on the Arduino Due board, for the purpose of collecting the probe readings. Uploaded to the board using the Arduino IDE desktop application. 
- [`ProbeOperation.ino`](<./Probe Drivers/ProbeOperation/ProbeOperation.ino>): For normal probe operation with functioning potentiometer
- [`TestPotentiometer.ino`](<./Probe Drivers/TestPotentiometer/TestPotentiometer.ino>): Prints outputs every second for general diagnostics
- [`ManualProbeOperation.ino`](<./Probe Drivers/ManualProbeOperation/ManualProbeOperation.ino>): Returns outputs for each button press, with manually recorded depths. The current version for running experiments (SHE, May 2023). 

To record the data straight to PC (as these drivers anticipate), the probe will be operated via PuTTY (or similar), rather than the Arduino serial command. Within the PuTTY settings, "Logging" will be on to a .csv file.

### MATLAB Processing
- [`probe_read.m`](./probe_read.m) : Matlab function for reading in probe data
- [`ProbeParse.m`](./ProbeParse.m) : Matlab script placed inside folder containing each probe's data which parses data to probe_read.m. Outputs a _.mat_ file containing a structure relateing to the density profile; plots the profiles; and prints the fitting parameters
- [`calc_N2.m`](./calc_N2.m) : Calculate the $N^2$ buoyancy frequency profile for a given density profile


## User Guide
### Installation
With the Arduino board plugged in via USB cable, and powered on open the .ino code in Arduino IDE. Within Tools>Board, find the relevant board, likely an Arduino Due, which may require the board manager to find. Then select the Port (in Tools>Port) corresponding to the Arduino board, it'll be a COM followed by a number, which you'll need later, and which should stay the same for this computer, unless major changes are made. To install the .ino driver, run the Sketch>Upload, this will initially check the code, and then upload it to the Arduino board. In the bottom of the screen, the output window should log progress and/or failure. Once this is successfully completed, the code is uploaded to the arduino, and should stay there until you need to update it. NOTE, this isn't in practice true, for some reason it periodically needs re-uploading, ensuring the Arduino is powered down in between experiments (power cable off, USB unplugged), seems to help this.  

### Logging the data
Probe operation differs a little between the probe driver files, and users should follow prompts in the serial command output to identify these differences. Usually, inputting R starts recordings, S stops them, P prints the output to the serial output, D deletes the most recent save to the SD card. The physical button usually initialises the readings, and the user should wait for a response from the serial command before commencing. In manual operation cases, the button is also used to take readings. 

Plug in the Arduino board, and set the probes at the surface of the water.
1) Open PuTTY, and you will be presented initially with the Configuration screen.
2) In this, you'll want to set `Connection type` to `Serial`, and set the Serial line to the COM port number previously identified for this board. 
3) Assuming you want to save the output to the computer, open the "Logging" tab, and set Session logging to "Printable output", Log file name to: "D:\path-to-probe-data\&D&M&Y\Probes_&T.csv". I also set "What to do if the log file already exists" to "Always overwrite it". 
4) Return to the Session tab, and enter a meaningful name "PROBES" into the Saved Sessions box, and "Save". Now you can simply click that setting each time and load it in future sessions. 
5) Click "Open" to start the session"
6) Carry out data collection as indicated onscreen/in the driver-specific documentation
7) Once completed, ensure the data has been printed through the Putty screen, either live, or by pressing P to print after the data collection. Now you can close the window, and look at the data in MATLAB. 




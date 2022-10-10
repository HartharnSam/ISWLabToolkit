% Routine for comparing the initial density profiles of each probe.
% Uses hydrometer readings to calibrate the probes
clearvars; close all; clc;
addpath('C:\Users\b5006861\OneDrive - Newcastle University\02_PhD_Project\01_PhD_Organisation\04_LaboratorySetup\01_SourceCodes\Current\probe_read');

%% user inputs
input_filename = 'D10.TXT';
final = true; % If confirmed final probe reading set to true

travel = 0.146;         % Probe travel distance
total_depth = 0.30;    % total depth
density_lims = [1029 1046];     % Density in [top bottom] layer from hydrometer
cut = [0 0];            % number of points to cut from [top bottom] of reading
% fix the percentage used for fit below!!

% Declare channels used   !!!!   be sure to do this   !!!!
    % pot-meter is on column 1 (channel 0)
    % so every other channels have a value of 1 higher
probes = [2];   % list of probe channels (ie. probes = [2 3 4];)
depth_shift = 0.00;     % Depth shift if dropped after this reading
pot_top_height = 0.30; % Probe height at top


%% Call function
[raw_data, fitted_data] = probe_read(input_filename, travel, total_depth, density_lims, cut, probes, depth_shift, pot_top_height);


%% Save data
if final == true
    calibr = struct('travel', travel, 'TotalDepth', total_depth, 'probes', probes, 'PotentiometerRestHeight', pot_top_height, 'HydrometerDensityLimits', density_lims, 'DataCut', cut, 'originalFile', input_filename);
    data = struct('RawData', raw_data, 'FittedData', fitted_data, 'Calibration', calibr);
    [~, dir_date, ~] = fileparts(pwd);
    fnm = [dir_date, '.mat'];
    save(fnm, 'data');
    
    print(gcf, [dir_date, '.png'], '-dpng')
end


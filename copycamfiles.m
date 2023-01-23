function copycamfiles(dateString, run_number, isNewWCS, isReprocess)
%COPYCAMFILES - Copies camera files from the original directory (D:) to
%OneDrive, and does processing of the ptv files, if needed copies new
%wcs-related status files into a shared directory
%
% Inputs:
%    dateString - 'Today', or the current date (e.g. '190122')
%    run_number - Number of the run (default 1)
%    isNewWCS - Boolean if needs to update WCS
%    isReprocess - Re-runs particle_tracks
%
% Other m-files required: particle_tracks
% Subfunctions: none
% MAT-files required: none
%
% See also: setup.sh,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-Jan-2022; Last revision: 19-Jan-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

if strcmpi(dateString, 'TODAY')
    dateString = datestr(now, 'ddmmyy');
end
if nargin < 2
    run_number = 1;
end
if nargin < 3
    isNewWCS = false;
end
if nargin < 4
    isReprocess = true;
    
end
digiflowstartup;
%% Fixed parameters
% Identify the current PC/Cameras
listCompNames = {'19S-STB-48040','18S-STB-89472',  '18S-STB-48039'};
list_CamNames = {'CamA', 'CamB', 'CamC'};
curr_ComputerName = getenv('COMPUTERNAME');
curr_CompNumber = find(strcmp(curr_ComputerName, listCompNames));
curr_CamName = list_CamNames{curr_CompNumber};

% Set the Directories
D_CamDir = 'D:\Sam\IceCovWaters\Camera_Data\';
OD_CamDir = 'C:\Users\b5006861\OneDrive - Newcastle University\02_PhD_Project\04_Ice_Covered_Waters\02_Raw_data\01_CameraData\';

% Change folder if the second run of day
if run_number == 2
    D_CamDir = [D_CamDir, '\2_'];
    OD_CamDir = [OD_CamDir, '\2_'];
end

%% Collate filenames
new_cam_dir = fullfile(OD_CamDir, dateString, curr_CamName);
orig_cam_dir = fullfile(D_CamDir, dateString);

%% Copy Output Data Files
copyfile([orig_cam_dir, '\output_*.dfi'], new_cam_dir);
try
    outputs = ls([orig_cam_dir, '\output_*.dfi']);
    im = dfireadvel([orig_cam_dir, '\', outputs(1, :)]);
catch
    warning('Output image may not be uncompacted/uncompressed')
end

copyfile([orig_cam_dir, '\DigiFlow_Status.dfs'], new_cam_dir);
copyfile([orig_cam_dir, '\wcs_first_cut.dfi'], new_cam_dir);
copyfile([orig_cam_dir, '\first_cut.dfm'], new_cam_dir);

try
    copyfile([orig_cam_dir, '\ts_*.dfi'], new_cam_dir);
    try
        outputs = ls([orig_cam_dir, '\ts_*.dfi']);
        im = dfireadvel([orig_cam_dir,'\', outputs(1, :)]);
    catch
        warning('Timeseries image may not be uncompacted/uncompressed')
    end
end
disp(dateString)
disp('Output Files Copied')

if isNewWCS % Copy over new wcs to the parent (Holding) folder
    copyfile([orig_cam_dir, '\DigiFlow_Status.dfs'], [D_CamDir, '\Settings\']);
    copyfile([orig_cam_dir, '\CoordinateSystems.log'],  [D_CamDir, '\Settings\']);
    copyfile([orig_cam_dir, '\DigiFlow_Dialogs.dfs'], [D_CamDir, '\Settings\']);
    copyfile([orig_cam_dir, '\wcs_first_cut.dfi'],  [D_CamDir, '\Settings\']);
    disp('WCS Files Updated')
end

%% Copy PTV data if required
        copyfile([orig_cam_dir, '\piv_hr_*.dfi'], new_cam_dir);

try
            outputs = ls([orig_cam_dir, '\piv_hr_*.dfi']);
            im = dfireadvel([orig_cam_dir,'\', outputs(1, :)]);
        catch
            warning('PIV image may not be uncompacted/uncompressed')
        end
if curr_CompNumber == 3
    % Copy raw DigiFlow (.txt) output
%    copyfile([orig_cam_dir, '\ptv_*.txt'], new_cam_dir);
%    copyfile([orig_cam_dir, '\particles.dfd'], new_cam_dir);
    
    % Then do the .mat file output too
    if isReprocess || (exist([new_cam_dir, '/ptv_tracks.mat'])~=2)
     %   particle_tracks(orig_cam_dir);
    end
    try
        copyfile([orig_cam_dir, '\ptv_tracks.mat'], new_cam_dir);
        disp('PTV Files Copied')
    catch
%        copyfile([orig_cam_dir, '\ptv_*.dfi'], new_cam_dir);
        try
            outputs = ls([orig_cam_dir, '\ptv_*.dfi']);
            im = dfireadvel([orig_cam_dir,'\', outputs(1, :)]);
        catch
            warning('PIV TS image may not be uncompacted/uncompressed')
        end
        disp('PIV Files Copied')
    end
else
    % Copy over High res PIV
    try
        %copyfile([orig_cam_dir, '\piv_*.dfi'], new_cam_dir);
        warning('correct to piv_hr')
        try
            outputs = ls([orig_cam_dir, '\piv_hr_*.dfi']);
            im = dfireadvel([orig_cam_dir,'\', outputs(1, :)]);
        catch
            warning('PIV image may not be uncompacted/uncompressed')
        end
    catch
        warning('No PIV Images detected')
    end
    
    % Copy over PIV timeseries
    try
        copyfile([orig_cam_dir, '\piv_ts.dfi'], new_cam_dir);
        try
            outputs = ls([orig_cam_dir, '\piv_ts.dfi']);
            im = dfireadvel([orig_cam_dir,'\', outputs(1, :)]);
        catch
            warning('PIV TS image may not be uncompacted/uncompressed')
        end
    catch
        warning('No PIV TS Images detected');
    end
    
end
disp('Done')
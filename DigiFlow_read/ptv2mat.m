function ptv2mat(filename, max_t, isPlot)
%PTV2MAT - reads in ptv data from each timestep incrementally, and
%collates for each particle track, saves as ptv_tracks.mat
%Can also track numbers vs locations over time if isPlot is "true" (doesn't save this)
%
% Other m-files required: ptv_read, completion
% Subfunctions: none
% MAT-files required: none
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Nov-2020; Last revision: 08-Feb-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%% Run the function
if nargin < 3
    isPlot = false;
end

max_files = numel(dir([filename, '/ptv_bas*']))-1;
if nargin<2 || max_files<max_t
    max_t = max_files;
end

for i = 1:max_t
    output = ptv_read(filename, i-1, 'full');
    emptyoutput = [];
    for j = 1:length(output)
        if ~isempty(output{j})
            emptyoutput(end+1) = j;
        end
    end
    if i == 1
        
        % Put something into the ptv structure to identify variables
        ptv.Variables = fieldnames(output{emptyoutput(1)});
        ptv.Variables{end+1} = 'x_Area';
        ptv.Variables{end+1} = 'y_Area';
        ptv.data = cell(1, 1);
    end
    
    for j = 1:length(emptyoutput)
        ind = emptyoutput(j);
        try
            if (output{ind, 1}.idTrack > length(ptv.data)) || isempty(ptv.data{output{ind, 1}.idTrack})
                ptv.data{output{ind, 1}.idTrack} = nan(max_t, length(ptv.Variables));
            end
            
            % Calculate the area centroids
            output{ind}.x_Area = output{ind}.x - output{ind}.mismatchx;
            output{ind}.y_Area = output{ind}.y - output{ind}.mismatchy;
            
            % Store each particle's data
            
            ptv.data{output{ind, 1}.idTrack}(i,:) = cell2mat(struct2cell(output{ind}));
            
        catch
           % error(['Error on processing timestep ' num2str(i)])
        end
        if isPlot
            text(output{ind, 1}.x, output{ind, 1}.y, num2str(output{ind, 1}.idTrack));
            hold on
        end
        
    end
    
    if isPlot
        drawnow
        hold off
        xlim([0 1000])
        ylim([0 1000])
        pause(.1)
        clf
    end
    completion(i, max_t);
end
ptv.n_particles = length(ptv.data);
ptv.n_timesteps = i;

save(fullfile(filename, 'ptv_tracks.mat'), 'ptv');


function ptv = ptv_read(filename, framenumber, importtype)
%PTV_READ - Imports PTV data, produced by the ptv_export.dfc code
%
%
% Inputs:
%    filename - Description
%    framenumber - Description
%    importtype - Description
%
% Outputs:
%    ptv - Description
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 13-Apr-2021; Last revision: 13-Apr-2021
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Import Basic Data
% Set up filename
basic_filepart = strrep('ptv_basic_####.txt', '####', sprintf('%04d', framenumber));
basic_fname = fullfile(filename, basic_filepart);

% Set table import options
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["x", "y", "u", "v", "ID"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "skip";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "x", "TrimNonNumeric", true);
opts = setvaropts(opts, "x", "ThousandsSeparator", ",");

% Import the data
ptvbasic = readtable(basic_fname, opts);
ptvbasic.ID = ptvbasic.ID + 1; % Add 1 to avoid 0s

ptvshort = cell(max(ptvbasic.ID), 1);
for i = 1:length(ptvbasic.ID)
    ii = ptvbasic.ID(i);
    ptvshort{ii, 1}.x = ptvbasic.x(i);
    ptvshort{ii, 1}.y = ptvbasic.y(i);
    ptvshort{ii, 1}.u = ptvbasic.u(i);
    ptvshort{ii, 1}.v = ptvbasic.v(i);
end
ptv = ptvshort;
%clear ptvbasic opts

%% Import complex data
if nargin < 3 || strcmpi(importtype, 'full') % Conditions to import full data
    clear opts;
    full_filepart = strrep('ptv_full_####.txt', '####', sprintf('%04d', framenumber));
    full_fname = fullfile(filename, full_filepart);
    opts = delimitedTextImportOptions("NumVariables", 29);
    
    opts.DataLines = [2, Inf];
    opts.Delimiter = " ";
    % Specify file level properties
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "skip";
    
    opts.VariableNamesLine = 1;
    % Specify column names and types
    opts.VariableNames = ["startFrame", "startIndex", "endFrame", "endIndex", "x", "y", "iTo", "iFrom", "idTrack", "volume", "vXYCorrelation", "vRMSSizex", "vRMSSizeyy", "area", "aXYCorrelation", "aRMSSizex", "aRMSSizey", "boxxMin", "boxxMax", "boxyMin", "boxyMax", "predictNextx", "predictNexty", "mismatchx", "mismatchy", "nEdgePoints", "sleeping", "threshold", "reject"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Import the data
    ptvfull = readtable(full_fname, opts);
    ptvfull.idTrack = ptvfull.idTrack +1;
    
    for i = 1:height(ptvbasic)
        trackID = ptvbasic{i, 5};
        ptv_loc = ptvfull.idTrack == trackID;
        
        for j = 1:size(ptvfull, 2)
            if ptvfull.x(i) == 0 && ptvfull.y(i) == 0 % Checks for "fake particles"
                ptv{trackID, 1}.(ptvfull.Properties.VariableNames{j}) = NaN;
            else
                ptv{trackID, 1}.(ptvfull.Properties.VariableNames{j}) = ptvfull{ptv_loc, j};
            end
            ptv{trackID, 1}.idTrack = trackID;
        end
    end
    %ptv = ptvlong;
end
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
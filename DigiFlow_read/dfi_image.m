clc; clearvars; close all;
% dfi_image - Plots dfi outputs from DigiFlow, by reading in DigiFlow Image
% files and output to an avi movie or single frame image
% Note: input .dfi must be uncompressed and uncompacted

% other m-files required: subaxis.m v1.1.0.0
%                         (uk.mathworks.com/matlabcentral/filexchange/3696-subaxis-subplot)
%                        plot_dfi

%% DEVELOPMENT
% Currently loops mixed up, so won't be able to do a video with multiple
% frames

%% User Set Variables
% Input Files
InputDirectory = pwd ;
fnm(1) = "piv_####";          % Filename of generic PIV .dfi image
fnm(2) = "output_####";
m(1) = 597;                  % First frame number
m(2) = 597;                  % First frame number for second panel
%m(3) = 450;
n = 0100;                     % Number of frames

% Output Files
outputDirectory = pwd;
outputFile = "dissipation";          % Output video filename
outputFileType = 'mp4';                % Image save type

BackgroundType(1) = "u";
BackgroundType(2) = "image";
ForegroundType(1) = "streamline";
ForegroundType(2) = "none";

dn = 52;                              % mesh size to display data
%                                  Recommended:
%                                   - Streamline - 52
%                                   - VelArrows - 12
vertical_stack = true;

%% Rarely changed variables
fnm_placeholder = "####";
q_scale = .20;                  % size of quiver arrows in plot for unit magnitude vector

% Slope in video
x1 = [6.36+(-6.9308e-04 *1023) 6.36]-.5;                 % Min and Max World x Coordinates of slope
y1 = [0.05 .1];
n_panels = length(m);

% Settings for stacked plots
subaxis_settings = ("'SpacingVert', .02, 'Margin', .09");

ArgsIn = struct('colormap', singlecycle, 'q_scale', q_scale,...
    'Resolution', dn); %, 'x1', x1, 'y1', y1);

%% Compute Variables

if vertical_stack
    ni = n_panels; mi = 1;
    x_axis_labelled = n_panels; % Subplot number of axis for which time label will be labelled
    y_axis_labelled = 1:n_panels;
else
    ni = 1; mi = n_panels; %#ok<UNRCH> 
    y_axis_labelled = 1;
    x_axis_labelled = 1:n_panels;
end

fig = figure;
if strcmp(outputFileType, 'avi')
    vid = VideoWriter(outputFile); % Output file
elseif strcmp(outputFileType, 'mp4')
    vid = VideoWriter(outputFile, 'MPEG-4'); % Output file
end

%% Read in file
Grids = cell(n_panels, 1);
for t = 1:n
    for ii = 1:n_panels
        fnm_0 = fullfile(InputDirectory, join([fnm(ii), '.dfi'], ''));
        
        % Sort out the axes
        i = ceil(ii/mi);
        j = mod(ii-1, mi)+1;
        
        eval(join(['ax(ii) = subaxis(ni, mi, j, i, ', subaxis_settings, ');'], ""));
        
        if any(ii == x_axis_labelled)
            is_xaxisLabelled = true;
        else
            is_xaxisLabelled = false;
        end
        
        if any(ii == y_axis_labelled)
            is_yaxisLabelled = true;
        else
            is_yaxisLabelled = false;
        end
        
        fnm_in = strrep(fnm_0, fnm_placeholder, sprintf('%04d',m(ii)+t-1));
        if t == 1
            [~, ~, Grids{ii}] = plot_dfi(fnm_in, BackgroundType(ii), ForegroundType(ii), ArgsIn, [], ax(ii));
        else
            plot_dfi(fnm_in, BackgroundType(ii), ForegroundType(ii), ArgsIn, Grids{ii}, ax(ii));
        end
        
        if n>1
            if t == 1 && ii == 1
                vid.FrameRate = Grids{ii}{1}.FrameRate;
                open(vid);
            end
            frame = getframe(fig);
            writeVideo(vid,frame);
        end
        
        
    end
end

if n > 1
    close(vid);
else
    print(outputFile,strcat(['-d', outputFileType])) % saves image as eps - can comment out if not needed
end


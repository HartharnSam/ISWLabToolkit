function combine_digiflow_avi(run_name)
%COMBINE_DIGIFLOW_AVI - Combines DigiFlow .avi movies from multiple cameras
%to a single .mp4 file. 
%Matches timings, location and attempts to match brightness overall
%
% Syntax:  digiflow_avi_combine(run_name)
%
% Inputs:
%    run_name - Name of experimental run (e.g. 090120)
%
%
% Other m-files required: dfireadvel, singlecycle (DigiFlow colormap)
% Subfunctions: none
% MAT-files required: none
%
% See also: %COMBINE_DFI_IMAGES
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@ncl.ac.uk
% August 2020; Last revision: 12-Aug-2020
% MATLAB Version: 9.6.0.1114505 Release: 2019a Service Pack: Update 2

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Get into the correct directory
% Will need to set path to camera data (stored in a directory with
% sub-directories for each camera named Camera1, Camera2 etc.)
cd(run_name)
cameras_list = dir('Camera*'); % Identifies the cameras available
n_cameras = length(cameras_list);

%% Set up coordinate systems for each image

xs = nan(n_cameras, 2000); ys = xs; x2s = xs; % Preallocate arrays
x0 = nan(1, n_cameras); y0 = x0; dx = x0; dy = x0; nx = x0; ny = x0;
x = nan(n_cameras, 2); y = x;

for i = 1:n_cameras
    piv_list = dir([cameras_list(i).name,'\piv_*']); % Loads a piv file, as
    % these contain all coordinate data for dfireadvel, may want to be any generic .dfi image?
    
    im=dfireadvel([piv_list(1).folder, '\', piv_list(1).name]); %File name, read in uncompressed AND uncompacted piv file
    [x(i, :), y(i, :), dx(i), dy(i), nx(i), ny(i), x0(i), y0(i)]  = dfi_grid_read(im);
    [~, ~, ~, ~, ~, ~, ~, ~, xs(i, 1:nx(i)), ys(i, 1:ny(i))] = dfi_grid_read(im);
end

%% Match times (time as number of frames)
start_time = importdata('limits.txt')'; % Start time of each clip relative to a known 0 point (e.g. lights on)
start_shift = max(start_time) - start_time; % time from start of each clip to when combined clip starts

% Preallocations
v = cell(1, 2);
NumFrames = nan(1, n_cameras);
end_times = nan(n_cameras, 1);
for i = 1:n_cameras
    v{i} = VideoReader([cameras_list(i).folder, '\',cameras_list(i).name, '\', run_name, '.avi']); % Load clips to access stats
    if isfield(v{i}, 'NumFrames')
        NumFrames(i) = v{i}.NumFrames; % If possible load number of frames directly
    else
        NumFrames(i) = v{i}.Duration.*v{i}.FrameRate;  % if not, calculate it
    end
    
    end_times(i,1) = NumFrames(i)+start_time(i); % Absolute number of frames from 0 point to end of each clip
end
end_shift = min(end_times) - start_time; % time from start of each clip to when combined clip ends
length_combined = end_shift(1) - start_shift(1);

% Preallocate
sequence = nan(n_cameras, length_combined+1);
n2x = nan(1, n_cameras);

for i = 1:n_cameras
    sequence(i, :) = (start_shift(i):end_shift(i))+1; % identify the frames to be used for each camera
    
    if verLessThan('matlab', '8.9')
        
        frame_temp = read(v{1, i}, sequence(i, 1)); % load frame into temporary holder
        frameA = frame_temp(:, :, 1);
        n2x(i) = size(frameA, 2);
    else
        frame_temp = im2gray(read(v{1, i}, sequence(i, 1))); % load frame into temporary holder
        
        n2x(i, :) = size(frame_temp);
        R(i) = imref2d(n2x(i, :), -x(i, :), y(i, [2 1]));
    end
    
    x2s(i,1:n2x(i))= linspace(x0(i), x(i, 2), n2x(i));
    clear frameA
end


for i = 2:n_cameras
    x_overlapA(i-1, :) = find(x2s(i-1,:) < x0(i));
    x_overlapB(i-1, :) = find(x2s(i, :) > x(i-1, 2));
end
%% Match intensities

intensity_A = nan(n_cameras-1);
intensity_B = nan(n_cameras-1);
for i = 2:n_cameras
    if verLessThan('matlab', '8')
        
        frame_temp = read(v{1, i}, sequence(i, 5));
        frameB(:, :) = frame_temp(:, :, 1);
        frame_temp = read(v{1, i-1}, sequence(i-1, 5));
        frameA(:, :) = frame_temp(:, :, 1);
    else
        frameB = im2gray(read(v{1, i}, sequence(i, 5)));
        frameA = im2gray(read(v{1, i-1}, sequence(i-1, 5)));
        
    end
    intensity_A(i-1) = sum(frameA(:, x_overlapA(i-1, :)), 'all')/((x2s(i,x_overlapA(1)) - x(i, 2)).*dy(i)*ny(i));
    intensity_B(i-1) = sum(frameB(:, x_overlapB(i-1, :)), 'all')/((x2s(i,x_overlapA(1)) - x(i, 2)).*dy(i)*ny(i)) ;
    
end
intensity_ratio = intensity_A./intensity_B;
intensity_ratio = [1 intensity_ratio];
%% Do for whole movie
vi = VideoWriter('Combined.mp4', 'MPEG-4');

open(vi);
for i = 1:length_combined
    if verLessThan('matlab', '8')
        subaxis(1, 1, 1, 'Margin', 0);
    
    for ii = 1:n_cameras
        xlim(-[max(max(x)) min(min(x))])
        frame_ind = sequence(ii, i);
        frame = read(v{1, ii}, frame_ind)*intensity_ratio(ii);
        imagesc((-x(ii, :)), flip(y(ii, :)), (frame(:, :, 1)));
        colormap(singlecycle)
        daspect([1 1 1])
        
        hold on; clear frame
    end
    hold off
    set(gca, 'nextplot', 'replacechildren');
    f = getframe(gcf);
    else
                frame_ind = sequence(1, i);
        frame = im2gray(read(v{1, 1}, frame_ind))*intensity_ratio(1);
        frame_RGBA = ind2rgb(frame, singlecycle);
        
        frame_ind = sequence(2, i);
        frame = im2gray(read(v{1, 2}, frame_ind))*intensity_ratio(2);
        frame_RGBB = ind2rgb(frame, singlecycle);
    C = imfuse(frame_RGBA, R(1), frame_RGBB, R(2), ...
        'blend', 'Scaling', 'none');
    end
    
    writeVideo(vi, f);
    
    disp([num2str(i/length_combined *100), '% Complete'])
    
end
close(vi);

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
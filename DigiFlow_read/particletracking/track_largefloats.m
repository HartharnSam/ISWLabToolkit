%TRACK_LARGEFLOATS - Taking an ice_ts.mat Hovmoller plot, identifies
%the edges of the float (NaN Regions) and with user input, identifies the
%float location through time. 
%
% Other m-files required: FiniteDiff, smooth
% Subfunctions: none
% MAT-files required: ICE_TS
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% Oct-2022; Last revision: 23-Jan-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clc; clearvars; close all;

% User parameters
float_width = 1.2; % Length of the float
c_isw = 0.1; % Indicative speed of the wave, does not need to be accurate
n_floats = 1; % Number of float paths the code will record. NB, if the float passes off screen and returns, it will count as 2. 

cutoff_percentage = .02; % Percentage of edge of frame to cut off

%% Load in and plot initial data
% Load
load('ice_ts.mat');
data = im.cdata(:, :, 1);
% Plot
figure(1);
imagesc(grid.x, grid.y, data);
xlabel('x'); ylabel('t')
colormap('gray')
title('Initial Data')

%% Detect float edges
% Pre-allocate arrays
xpeaks = nan(1, grid.ny);
trans_points = [];
actual_track = nan(grid.ny, n_floats);

x_edges = [grid.xi(round(grid.nx*cutoff_percentage)), grid.xi(round(grid.nx*(1-cutoff_percentage)))];
% For each frame, do edge detection & validation
for ii = 1:grid.ny
    % uncomment for blow by blow account of the algorithm
    %figure(2);
    %plot(grid.xi, data(ii, :));
    %title('Frame-by-Frame Account')

    %Identify all the edges by finding "NaN Crossing points"
    nan_inds = find(isnan(data(ii, :)));
    nnan_inds = find(~isnan(data(ii, :)));

    edge_inds = nan_inds(diff([0 nan_inds])>1); % Put in places where it moves from not-NaN to NaN
    edge_inds = [edge_inds nnan_inds(diff([0 nnan_inds])>1)-1]; % And the places where it moves from NaN to not-NaN
    xpeaks_tmp = grid.xi(sort(edge_inds)); % Get corresponding x coordinates
    xpeaks_tmp(xpeaks_tmp>max(x_edges) | xpeaks_tmp < min(x_edges)) = NaN; % NaN if at the edge of frame

    %% Path splitting
    % Identify where lines are non-continuous & make a note
    if ii > 2 && isnan(xpeaks(ii-1)) && isnan(xpeaks(ii-2))
        trans_points = [trans_points ii];
    end

    % Identify if paths are not the same path & split
    if ii > 1 && ~isnan(xpeaks(ii-1)) && length(xpeaks_tmp)>=1 % check if there's adjoining points (t-1)
        region = c_isw*grid.dy; % criteria based on ISW speed
        xpeak_min = xpeaks(ii-1)-region;
        xpeak_max = xpeaks(ii-1)+region;
        xpeakind = find(xpeaks_tmp>xpeak_min & xpeaks_tmp<xpeak_max); % Identify which peak most closely matches the predicted path
        if length(xpeakind)>=1
            xpeaks_tmp(1) = xpeaks_tmp(xpeakind(1)); % Set the path location to the closest peak
        end

        if (xpeaks_tmp(1) < xpeak_min || xpeaks_tmp(1) > xpeak_max) || isnan(xpeaks_tmp(1)) % if the current location is within c_isw*dt
            % No closest path available, so identify a transition point
            warning(['On ii = ', num2str(ii), ' Switch Paths!'])
            xpeaks_tmp = NaN;
            trans_points = [trans_points ii];
        end
    end
    % Save this timestep's float location
    xpeaks_tmp = xpeaks_tmp(~isnan(xpeaks_tmp));
    if ~isempty(xpeaks_tmp)
        xpeaks(ii) = xpeaks_tmp(1);
    else
        xpeaks(ii) = NaN;
    end
end

% Make a plot of detected paths
figure(1);
hold on;
plot(xpeaks, grid.yi, 'b-');
legend('Initial paths')

%% Combine Paths and Identify centre point of each float
track_index = 1;
all_trans_points = [1 trans_points length(xpeaks)];
new_trans_points = [];
track_count = zeros(n_floats, 1);

% I think this gets rid of points that relate to really short paths
for ii = 1:length(trans_points)+1
    current_indexes = all_trans_points(ii):all_trans_points(ii+1);
    if abs(diff(current_indexes([1 end]))) > 3
        new_trans_points = [new_trans_points all_trans_points(ii)];
    else
        xpeaks(current_indexes) = NaN;
    end
end

new_trans_points = [new_trans_points length(xpeaks)];

for ii = 1:length(new_trans_points)-1
    current_indexes = new_trans_points(ii):new_trans_points(ii+1); % Identify the indices relating to the current path

    % Plot the current and the next path
    figure(3);
    imagesc(grid.x, grid.y, data);
    colormap('gray')
    hold on
    plot(xpeaks(current_indexes), grid.yi(current_indexes), 'r-');

    title('Manual Track Editing');
    ylabel('t (s)'); xlabel('x (m)')

    if ~(ii == length(new_trans_points)-1)
        next_indexes = new_trans_points(ii+1):new_trans_points(ii+2); % Identify the indices relating to the next path
        plot(xpeaks(next_indexes),grid.yi(next_indexes), 'b-');
    end

    hold off
    legend('Current Track', 'Later Track', 'location', 'bestoutside')
    drawnow;

    % Join tracks back together with user input
    track_index = input('Input Track Number: '); % Name each coherent track
    if ~isempty(track_index)
        disp('Indicate track connection - h (higher x location), l (lower x location)'); % Identify if the track is the front or rear of the float

        switch input('Input track direction: ', "s")
            case 'h'
                track_sign = -1;
            case 'l'
                track_sign = 1;
        end
        track_count(track_index) = track_count(track_index)+1;
        actual_track(current_indexes, track_index) = xpeaks(current_indexes)+(track_sign*float_width/2);
        nnan_ind = find(~isnan(actual_track(current_indexes, track_index)));
        if track_count(track_index) == 1
            track_correction(track_index) = actual_track(current_indexes(nnan_ind(end)), track_index);
        else
            correction = track_correction(track_index)-actual_track(current_indexes(nnan_ind(1)), track_index); % This correction joins tracks smoothly
            if abs(correction) < 0.02
                actual_track(current_indexes, track_index) = actual_track(current_indexes, track_index)+correction;
            end
            track_correction(track_index) = actual_track(current_indexes(nnan_ind(end)), track_index);
        end
    end
end


%% Smooth the actual track to remove jumps related to grid discretisation
u = actual_track*NaN;
u_sm = u;
for ii = 1:size(actual_track, 2)
    nnan_inds = find(~isnan(actual_track(:, ii)));
    actual_track(nnan_inds(1):nnan_inds(end), ii) = smooth(grid.yi(nnan_inds(1):nnan_inds(end)), ...
        actual_track(nnan_inds(1):nnan_inds(end), ii), .05);
    %% Calculate the float velocity
    nnan_inds = find(~isnan(actual_track(:, ii)));
    Dmat = FiniteDiff(grid.yi(nnan_inds), 1, 2, false, ~(max(diff(nnan_inds))>1));
    u(nnan_inds, ii) = Dmat*actual_track(nnan_inds, ii);
    start_end_ind = nnan_inds(1):nnan_inds(end);
    u_sm(start_end_ind, ii) = smooth( grid.yi(start_end_ind), u(start_end_ind, ii), .05,  'rlowess');

end

%% Plot the final tracks
figure(4);
imagesc(grid.x, grid.y, data);
colormap('gray')
hold on
plot(actual_track, grid.yi, 'b-x')
xlabel('x (m)'); ylabel('t (s)')
title('Finalised Particle Tracks')

figure(5);
plot(grid.yi, u, 'k-');

hold on
plot(grid.yi, u_sm, 'r-');
title('Particle Velocities')
xlabel('t (s)'); ylabel('u (m/s)');

%% Save the structure
for ii = 1:n_floats
    ptv.data{ii}(:, 1) = actual_track(:, ii);
    ptv.data{ii}(:, 3) = u(:, ii);

end

ptv.n_timesteps = grid.ny;
ptv.n_particles = n_floats;
ptv.Variables = {'x', 'y', 'u', 'v', 'startFrame', 'idTrack', 'startIndex', 'endFrame', 'endIndex', 'iTo', 'iFrom', 'area'};
save('ptv_tracks_compiled.mat', 'ptv');

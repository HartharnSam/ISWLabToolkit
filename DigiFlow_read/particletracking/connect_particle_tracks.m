%CONNECT_PARTICLE_TRACKS - Links together particle tracks created by
%ptv2mat, or track_largefloats, using user input. Then re-calculates
%velocities
%
% Other m-files required: FiniteDiff, smooth (Curve Fitting Toolbox),
% dfi_grid_read, dfireadvel
% MAT-files required: ptv_tracks_compiled
%
% See also: CALC_FLOAT_HOVMOLLER
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% Oct-2022; Last revision: 23-Jan-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

clearvars; close all; clc;

%% Bring together multiple ptv files
CamA = load('../CamA/ptv_tracks_compiled.mat');
CamB = load('../CamB/ptv_tracks_compiled.mat');
CamC = load('ptv_tracks_compiled.mat');

n_times = max([CamA.ptv.n_timesteps CamB.ptv.n_timesteps CamC.ptv.n_timesteps]);
%n_times = max([CamA.ptv.n_timesteps CamB.ptv.n_timesteps]);
%n_times = max([CamA.ptv.n_timesteps]);

new_ptv.data = CamA.ptv.data;
%CamA.ptv.n_particles = 3;
count = CamA.ptv.n_particles;

new_ptv.data(count+[1:CamB.ptv.n_particles]) = CamB.ptv.data;
count = count+CamB.ptv.n_particles;

CamC.ptv.n_particles = 2;
new_ptv.data(count+[1:CamC.ptv.n_particles]) = CamC.ptv.data;
count = count+CamC.ptv.n_particles;

new_ptv.n_particles = count;
new_ptv.n_timesteps = n_times;
new_ptv.Variables = CamA.ptv.Variables;

for ii = 1:count
    current_length = length(new_ptv.data{ii});
    new_ptv.data{ii}(current_length:n_times, :) = NaN;
end
ptv = new_ptv;
save('new_ptv_tracks_compiled.mat', 'ptv');
clearvars; close all;

%%
% Plot the Initial Data

dirname = '.';
load(fullfile(dirname, 'new_ptv_tracks_compiled.mat'), 'ptv'); % Load in the data
%load(fullfile(dirname, 'ptv_tracks_compiled.mat'), 'ptv'); % Load in the data

% Make a conversion from pixels to WCS
im = dfireadvel(fullfile(dirname, 'output_0000.dfi'));
Grid  = dfi_grid_read(im);
timestep = 1/30;
n_particles = ptv.n_particles;

time = (0:ptv.n_timesteps-1)*timestep;
variable_index = find(strcmpi(ptv.Variables, 'x'));
continuing = true;
while continuing
    close all;
    for i = 1:n_particles
        if ~isempty(ptv.data{i})
            %ptv.data{i}(:, variable_index) = interp1([1 Grid.nx], (Grid.x), ptv.data{i}(:, variable_index));
            active_inds = find(~isnan(ptv.data{i}(:, 1)));
            plot(time, ptv.data{i}(:, 1), 'DisplayName', num2str(i))
            hold on
        end
    end
    legend('Location', 'best')
    set(gcf, 'Position', [1221 376 560 420]);
    drawnow;
    % Identifying tracks to connect
    connecting_lines = input('Lines to connect (empty if none): ');
    if ~isempty(connecting_lines)

        % Connect the tracks
        active_inds_2 = find(~isnan(ptv.data{connecting_lines(2)}(:, 1)));
        ptv.data{connecting_lines(1)}(active_inds_2, :) = ptv.data{connecting_lines(2)}(active_inds_2, :);
        ptv.data{connecting_lines(2)} = [];

        % Then fill missing
        temp_data = fillmissing(ptv.data{connecting_lines(1)}, 'pchip', 1, 'EndValues', 'none');
        ptv.data{connecting_lines(1)} = temp_data;


    end
    continuing = input('Add another line? (True/False)');
end

close all
for i = 1:n_particles
    if ~isempty(ptv.data{i})
        %ptv.data{i}(:, variable_index) = interp1([1 Grid.nx], (Grid.x), ptv.data{i}(:, variable_index));
        active_inds = find(~isnan(ptv.data{i}(:, 1)));
        plot(time, ptv.data{i}(:, 1), 'DisplayName', num2str(i))
        hold on
    end
end
particles = find(~cellfun(@isempty, ptv.data));
ptv.data = ptv.data(particles);

set(gcf, 'Position', [1221 376 560 420]);
drawnow;

if input('Is Ok to save (true/false)? ')
    save(fullfile(dirname, 'ptv_tracks_compiled.mat'), 'ptv');
end

% Re-calculate velocities:
recalc = true;
actual_track = ptv.data{1}(:, 1);
Grid.yi = (0:ptv.n_timesteps)/30;
if recalc
    u = actual_track*NaN;
    u_sm = u;
    actual_track = ptv.data{1}(:, 1);
    nnan_inds = find(~isnan(actual_track));
    actual_track(nnan_inds(1):nnan_inds(end)) = smooth(Grid.yi(nnan_inds(1):nnan_inds(end)), ...
        actual_track(nnan_inds(1):nnan_inds(end)), .05);

    %% Calculate the float velocity
    nnan_inds = find(~isnan(actual_track));
    Dmat = FiniteDiff(Grid.yi(nnan_inds), 1, 2, false, ~(max(diff(nnan_inds))>1));
    u(nnan_inds) = Dmat*actual_track(nnan_inds);
    start_end_ind = nnan_inds(1):nnan_inds(end);
    u_sm(start_end_ind) = smooth(Grid.yi(start_end_ind), u(start_end_ind), .05,  'rlowess');
    ptv.data{1}(:, 3) = u_sm;
    save(fullfile(dirname, 'ptv_tracks_compiled.mat'), 'ptv');

end

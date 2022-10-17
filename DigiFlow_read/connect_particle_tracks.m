clearvars; close all; clc;

% Plot the Initial Data

filename = '.\CamC';
load('./CamC/ptv_tracks.mat', 'ptv'); % Load in the data
% Make a conversion from pixels to WCS
im = dfireadvel(fullfile(filename, 'output_0000.dfi'));
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
    figure_print_format(gcf)
    legend('Location', 'best')
    set(gcf, 'Position', [1221 376 560 420]);
    drawnow;
    % Identifying tracks to connect
    connecting_lines = input('Lines to connect');
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
set(gcf, 'Position', [1221 376 560 420]);
figure_print_format(gcf)

drawnow;

if input('Is Ok?')
    save('CamC\ptv_tracks_compiled.mat', 'ptv');
end
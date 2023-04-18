% CALC_FLOAT_HOVMOLLER - Construct Hovmoller for lab data with (optional?) 
% cutting of the ice region, and saves as ice_ts.mat 
% Reads in raw images frame by frame, identifies the float (using
% find_boundaries) and takes a row (ts_level) to construct the hovmoller. 
% 
% Other m-files required: dfireadvel, dfi_grid_read, find_boundaries,
% completion
% Subfunctions: none
% MAT-files required: none
%
% Development note: There is no question that this script would be
% MUCH faster if it were written in DigiFlow. Plotting does slow it down,
% but I think the human oversight is fairly important here. The concern
% over transfer to DigiFlow code came from the use of a gaussian filter
% and other bits and pieces in MATLAB script
%
% See also: TRACK_LARGEFLOATS,  FIND_BOUNDARIES
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 23-Jan-2023; Last revision: 23-Jan-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

clearvars; close all; clc;
%% Set Variables
im_ice_threshold = 85; % Brightness threshold for the float edge
ice_thickness = 0.04; % Maximum depth that ice would be detected at (in m)
ts_level = 510; % Row of pixels to be samples for the hovmoller

diagnostic_plot = true; % plots (intermittently) the hovmoller and raw image with ice removed for refernece. n.b. it slows down processing a lot

%% Run the code
% Preallocate arrays by reading in an example image and the grid
fnm_in = strrep('output_####.dfi', '####', sprintf('%04d',0));
im = dfireadvel(fnm_in);
grid = dfi_grid_read(im);

nt = length(dir('output_*.dfi'));
nx = im.nx;
ts_data = NaN(nt, nx);

for ii = 1:nt
    fnm_in = strrep('output_####.dfi', '####', sprintf('%04d',ii-1));
    im = dfireadvel(fnm_in);
    data = im.cdata(:, :, 1);
    % Add in ice and pycnocline boundaries
    % Threshold for brightness attributed to "Ice"
    II = find_boundaries(fnm_in, 'ice', ice_thickness, im_ice_threshold);
    y_ice = interp1(II.x, II.y_ice, grid.X(1, :));
    y_ice(isnan(y_ice)) = max(grid.y);

    ice_mask = false(im.nx, im.ny);
    for i = 1:im.nx
        ice_mask(i, grid.Y(:, i) > y_ice(i)) = 1;
    end
    data(ice_mask') = NaN;
    ts_data(ii, :) = data(end-ts_level ,: );

    if diagnostic_plot && mod(ii, 80)==0
        figure(1)
        subplot(1, 2, 1)
        pcolor(grid.X, grid.Y, data);
        colormap('gray'); caxis([0 255]);
        hold on
        yline(grid.Y(end-ts_level, 1), 'r');
        hold off
        subplot(1, 2, 2)
        if ii >1
            pcolor(ts_data);
            colormap('gray')
            caxis([0 255])
        end
        drawnow;
    end
    completion(ii, nt)
end

%% Save the thing
grid.ny = nt;
grid.dy = im.tStep;
grid.y = im.tFirst + [0 (nt-1)*im.tStep];
grid.yi = grid.y(1):grid.dy:grid.y(2);
[grid.X, grid.Y] = meshgrid(grid.xi, grid.yi);

clear im
im.cdata(:, :, 1) = ts_data;
settings = struct('im_ice_threshold', im_ice_threshold, 'ice_thickness', ice_thickness, 'ts_level', ts_level);

save('ice_ts.mat', 'im', 'grid', 'settings');


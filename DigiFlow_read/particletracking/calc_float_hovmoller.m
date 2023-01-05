% Construct Hovmoller for lab data with (optional?) cutting of the ice
% region
% Development note: There is no question that this script would be
% MUCH faster if it were written in DigiFlow. Plotting does slow it down,
% but I think the human oversight is fairly important here. The concern
% over transfer to DigiFlow code came from the use of a gaussian filter
% and other bits and pieces in MATLAB script

digiflowstartup;
clearvars; close all; clc;
im_ice_threshold = 28;
ice_thickness = 0.04;
ts_level = 480;

diagnostic_plot = true; % Slows it down a lot

% Preallocate arrays
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
        pcolor(grid.X, grid.Y, data);
        colormap('gray'); caxis([0 255]);
        hold on
        yline(grid.Y(end-ts_level, 1), 'r');
        hold off
        figure(2)
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
save('ice_ts.mat', 'im', 'grid');


function [im, x, y] = combine_im(filename)
%COLLATE_IM - Returns image, and grid for multiple images stitched
%together. 
% It interpolates, so it's not very efficient, but needed when
%we're calculating things based on it ~0.5s slow for three frames.
%
% Inputs:
%    filename - Cell array of filenames to stitch, e.g.
%    {'CamA\piv_0001.dfi', 'CamB\piv_0001.dfi'};
%
% Outputs:
%    im - image structure, like that created by dfireadvel
%
% Other m-files required: dfireadvel, dfi_grid_read
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 20-Feb-2023; Last revision: 20-Feb-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

cutoff = 8; % Amount to cut off the edges, which appears to be an effect on PIV images

%% Read in each image
im = cell(1, length(filename)); grid = im;
for ii = 1:length(filename)
    im{ii} = dfireadvel(filename{ii});
    grid{ii} = dfi_grid_read(im{ii});
    im{ii}.cdata([1:cutoff end-cutoff:end], :, :) = NaN;
    im{ii}.cdata(:, [1:cutoff end-cutoff:end], :) = NaN;

    % Identify the edges
    if ii == 1
        xmin = min(grid{ii}.x);
        xmax = max(grid{ii}.x);
        ymin = min(grid{ii}.y);
        ymax = max(grid{ii}.y);
    else
        xmin = min(xmin, min(grid{ii}.x));
        xmax = max(xmax, max(grid{ii}.x));
        ymin = min(ymin, min(grid{ii}.y));
        ymax = max(ymax, max(grid{ii}.y));
    end
end
% Re-create grids using grid spacing from image 1
if grid{1}.dx > 0
    new_x = xmin:grid{1}.dx:xmax;
else
    new_x = xmax:grid{1}.dx:xmin;
end
new_y = ymin:grid{1}.dy:ymax;

n_params = size(im{1}.cdata, 3);

[newX, newY] = meshgrid(new_x, new_y);
im2 = NaN(size(newX, 1), size(newX, 2), n_params);

% Interpolate data from each image into the collated image
for jj = 1:n_params
    data = newX*NaN;
    for ii = 1:length(filename)
        tmp_data = interp2(grid{ii}.X, grid{ii}.Y, im{ii}.cdata(:, :, jj), newX, newY);
        data(~isnan(tmp_data)) = tmp_data(~isnan(tmp_data));
    end
    im2(:, :, jj) = data;
end
clear im
y = new_y';
x = new_x';

im.cdata = im2;
im.x = new_x'; im.y = new_y';
im.nx = size(im.x, 1); im.ny = size(im.y, 1);
im.xWorldPerPixel = im.x(2) - im.x(1);
im.yWorldPerPixel = im.y(2) - im.y(1);
im.xOriginWorld = im.x(1); im.yOriginWorld = im.y(1);

end
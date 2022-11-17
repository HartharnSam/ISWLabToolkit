% Plot to try and figure out float tracking for large floats
% use on 180322_A
%
clc; clearvars; close all; 

ice_threshold = 180/255;

im = dfireadvel('ice_ts.dfi');
grid = dfi_grid_read(im);

imagesc(grid.x, grid.y, im.cdata(:, :, 1));
colormap('gray')
[brightnesses, x_inds] = sort(im.cdata(:, :, 1), 2);

xs = grid.xi(x_inds);

%inds = (brightnesses(:, end)>ice_threshold) && abs(brightnesses(:, end)-brightnesses(:, ));
hold on
ranges = range(xs(:, end-5:end), 2);
plot(ranges+4.7, grid.yi, 'bx')
inds = ranges < 0.01;
plot(ranges(~inds)+4.7, grid.yi(~inds), 'rx')

hold on
plot(xs(inds, end), grid.yi(inds));

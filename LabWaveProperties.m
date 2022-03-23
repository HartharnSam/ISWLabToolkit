function LabWaveProperties(filename, y_ref, colormap)
%LabWaveProperties - Takes user mouse click input to measure streamline
%displacement of digiflow timeseries image. Saves all to png image
%
% Inputs:
%    filename - filename as "ts_c200" (potentially including path), without
%    extension
%    y_ref - [Optional] Reference y line to plot (usually pyc depth)
%    colormap - [Optional] Colormap to use
%
% Other m-files required: dfireadvel; plot_dfi; mginput
% Subfunctions: none
% MAT-files required: none
%
% See also:  
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 15-Feb-2022; Last revision: 15-Feb-2022
% MATLAB Version: 9.10.0.1851785 (R2021a) Update 6

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
close all;
input_fnm = join([filename, ".dfi"], '');

% Load and plot data
if nargin > 2
    ArgsIn.colormap = colormap;
else
    ArgsIn.colormap = 'greyscale';
end
plot_dfi(input_fnm, 'TimeSeries', 'none', ArgsIn); hold on; 
xlabel('t (s)'); ylabel('y (m)')
set(gca, 'XDir', 'normal')

% Plot on reference line 
if nargin > 1
yline(y_ref, 'w-', 'LineWidth', 2); yline(y_ref, 'k--', 'LineWidth', 1);
end

% 
im = dfireadvel(input_fnm);
contour(im.x, im.y, im.cdata, 'k-')

% Pick points & Plot picked points
[x, y] = mginput(2); % Pick
plt = plot(x, y, 'r-x'); % Plot
datatip(plt, x(1), y(1)); datatip(plt, x(2), y(2), 'location', 'southwest'); % Show

% Display results to Command Window
fprintf('Streamline Traced : %.3f \n', y(1));
fprintf('Amplitude : %.3f m \n', -diff(y));
fprintf('Time min : %.3f s \n', x(2));

% Save figure with the data
exportgraphics(gcf, join([filename, '.png'], ''));

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
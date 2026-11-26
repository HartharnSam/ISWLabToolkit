%EXAMPLE_PLOT_TIMESLICES - Worked example of constructing ColumnProperties
% and calling plot_timeslices to compare a SPINS model run against lab data.
%
% NOTE: every column in ColumnProperties must have the same number of
% entries in Times (n_times is taken from column 1 and reused for all
% columns) — the two columns below are kept in sync for this reason.
%
% Other m-files required: plot_timeslices
% Subfunctions: none
% MAT-files required: none
%
% See also: plot_timeslices
% 25-Nov-2026; Last revision: 25-Nov-2026
% MATLAB Version: 26.1.0.3251617 (R2026a) Update 2

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
clc; clearvars; close all;

%% Column 1: SPINS model output, density field
ColumnProperties{1}.DirectoryName = '../02_Raw_data/model_run_01';
ColumnProperties{1}.Type = 'model';
ColumnProperties{1}.Times = [20 40 60];     % output numbers to plot
ColumnProperties{1}.Field = 'rho';
ColumnProperties{1}.xlimits = [0 5];
ColumnProperties{1}.FGType = 'vorty';       % overlay vorticity contours
ColumnProperties{1}.contourlevels = [-1 -0.5 0.5 1];

%% Column 2: Lab (DigiFlow) PIV data, vertical velocity field
ColumnProperties{2}.DirectoryName = {'../01_Lab_data/run_01/cam1'};
ColumnProperties{2}.Type = 'lab';
ColumnProperties{2}.Times = [20 40 60];     % must match column 1's time count
ColumnProperties{2}.Field = 'w';
ColumnProperties{2}.xlimits = [0 5];

%% Call plot_timeslices
c_range = 2*[-1 1];        % colour axis limits applied to both columns
labelled_columns = 2;      % only the lab column gets a colourbar
ylims = [-0.3 0];          % shared y-axis limits

figure;
ax = plot_timeslices(ColumnProperties, c_range, labelled_columns, ylims);
% ax{1}, ax{2} - axes handles for each column, one per time, for further
% customisation (e.g. ax{1}(1) is column 1's first-time axes)

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
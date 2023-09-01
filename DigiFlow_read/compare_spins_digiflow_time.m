function compare_spins_digiflow_time(spins_dir, spins_ii, digiflow_dir, digiflow_ii)
%COMPARE_SPINS_DIGIFLOW_TIME - Compare times for SPINS and DigiFlow outputs
%Plots the given timestep for the model and lab to match times up
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 17-Aug-2023; Last revision: 17-Aug-2023
% MATLAB Version: 9.12.0.2170939 (R2022a) Update 6

%Test environment
%  digiflow_dir = './Lab/110320/CamA/output_####.dfi';
%  digiflow_ii =  10;
%  spins_dir = './Model/NovakTank/260620_22';
%  spins_ii  = 72;
%/test environment

spinsstartup; digiflowstartup;

% Load in and plot digiflow part
    ax1 = subplot(2, 1, 1);

if ischar(digiflow_dir) % checks if one or more than one camera
    fnm_in = strrep(digiflow_dir, '####', sprintf('%04d',digiflow_ii));
    im = dfireadvel(fnm_in);
    grids = dfi_grid_read(im);
    imagesc(grids.x, grids.y, im.cdata(:, :, 1));
    xlims = sort(grids.x);

end
set(ax1, 'YDir', 'normal')
%colormap(ax1, singlecycle)
colormap(ax1, cmocean('balance'));
clim([-1 1].*0.1);
hold on;
grid on;

% load in and plot SPINS part
orig_dir = cd;
cd(spins_dir)
[x, z] = spinsgrid2d;
area(ax1, x(:, 1)-.6, z(:, 1)+.3, 'LineStyle', '-', 'FaceColor', 'w');
xlim(xlims)

params = spins_params;
x = x-params.L_adj -0.3;
rho = spins_reader_new('u', spins_ii);
ax2 = subplot(2, 1, 2);
pcolor(x, z, rho); 
colormap(ax2, cmocean('balance'));
xlim(xlims)
clim([-1 1].*0.1);
grid(ax2, 'on'); set(ax2, 'Layer', 'top')
cd(orig_dir)
end
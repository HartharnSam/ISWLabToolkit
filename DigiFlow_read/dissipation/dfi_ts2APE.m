function APE = dfi_ts2APE(filename, t0, t1, delta_rho, c0)
%DFI_TS2APE - Calculates an approximate APE from a column timeseries image
%Using the method of Boegman et al. (2005). JFM, 531, 159-180. doi:10.1017/S0022112005003915
%Identifying the pycnocline displacement via the pycnocline detection
%method in FIND_BOUNDARIES (Peter Sutherland).
%
% APE = cg * delta_rho * int_t0^t1(eta^2) dt
%
% Inputs:
%    filename - image file name (e.g. 'ts_712.dfi')
%    t0 - Lower limit for time for integration
%    t1 - Upper extent of time for integration
%    delta_rho - rho_2-rho_1
%    c0 - Wave propagation speed
%
% Outputs:
%    APE - Calculated APE
%
% Other m-files required: dfi_grid_read, dfireadvel, find_boundaries
% Subfunctions: none
% MAT-files required: none
%
% See also: FIND_BOUNDARIES
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Jan-2023; Last revision: 05-Jan-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

close all; 
im = dfireadvel(filename);
grid = dfi_grid_read(im);

pcolor(grid.xi, grid.yi, im.cdata);

II = find_boundaries(filename, 'pycnocline', [], [], false);

eta = smooth(II.x, II.y_interface, .1, 'rloess');
hold on
plot(II.x, eta, '-k')
ref_y = eta(1);
%ref_y = 0.2356;
eta = eta - ref_y;

%% Cut off t0 - t1
inds = find((II.x >= t0) & (II.x < t1));
eta = eta(inds);
x = II.x(inds);

% Integrate eta^2
int_eta2 = trapz(x, eta.^2);

% Calculate APE
APE = c0*9.81*delta_rho*int_eta2;

end
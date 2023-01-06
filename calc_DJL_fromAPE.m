function [wavelength, DJL] = calc_DJL_fromAPE(h_1, h_2, h_pyc, APE, rho_1, rho_2)
%CALC_DJL_fromAPE - Calculates a DJL solution to match the lab wave, based off the
%measured wave APE.
%
% Inputs:
%    h1 - Upper layer depth (depth to pyc centre in tanh profile)
%    h2 - Lower Layer depth
%    h_pyc - Thickness of pycnocline
%    amp - Amplitude of wave to find the wavelength for
%    rho_1 - Upper layer density (kg/m^3)
%    rho_2 - Lower Layer Density (kg/m^3)
%
% Outputs:
%    wavelength - Computed DJL Wavelength
%    DJL - Structure containing a host of interesting outputs
%
% Other m-files required: DJLES package (https://github.com/mdunphy/DJLES)
% Subfunctions: none
% MAT-files required: none
%
% See also: case_large_ape,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 08-Mar-2022; Last revision: 08-Mar-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
m_path = mpath;
addpath([m_path, 'DJLES'])

%% Specify/Load parameters
A = APE; 
L  = 14.0;  % domain width (m)
H  = h_1+h_2;  % domain depth (m)

relax = 0.15; 

% Specify the general density profile which takes d_d as a second parameter
delrho = abs(rho_1 - rho_2)/rho_2; %
a_d = delrho/2;
z0_d = h_1 + (h_pyc/2); % depth of pycnocline
frho=@(z,d_d) 1-a_d*tanh((z+z0_d)/d_d);
frhoz=@(z,d_d) -(a_d/d_d)*sech((z+z0_d)/d_d).^2; % derivative in z

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%% Find the solution %%%
d_d = h_pyc;
NXList = [64 128 256 512 1024];
NZList = [32 64  128 128 256];
rho  = @(z) frho(z, d_d);
rhoz = @(z) frhoz(z, d_d);
verbose = 0;

start_time = clock;
for ddindex = 1:length(NXList)
    % Resolution for this wave
    NX = NXList(ddindex);
    NZ = NZList(ddindex);

    % Density profile for this wave
    % Iterate the DJL solution
    djles_refine_solution
end

% Reduce epsilon, iterate to convergence
epsilon = 1e-5;
djles_refine_solution;

% Raise resolution, iterate to convergence
NX = 2048; NZ = 1024;
djles_refine_solution;

end_time = clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Save loopy things
amp = -wave_ampl;
% Compute and plot the diagnostics
djles_diagnostics
djles_plot
wavelength = djles_wavelength(eta, L);

% Save some loopy things

DJL.WaveAPE = A;
DJL.WaveAmp = -wave_ampl;
DJL.WaveC = c;
DJL.density = density; 
DJL.eta = eta; DJL.etax = etax; DJL.etaz = etaz;
DJL.u = u; DJL.ux = ux; DJL.uz = uz;
DJL.vorticity = vorticity;
DJL.w = w; DJL.wx = wx; DJL.wz =wz; 
DJL.WaveWavelength = wavelength; 
DJL.x = xc; DJL.z = zc;
save('DJL', 'DJL', 'uwave', 'c', 'x', 'z', 'density', 'L', 'wavelength');
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------


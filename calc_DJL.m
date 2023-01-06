function [wavelength, DJL] = calc_DJL(h_1, h_2, h_pyc, wave_amp, rho_1, rho_2)
%CALC_DJL - Calculates a DJL solution to match the lab wave, based off the
%measured wave amplitude.
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
% See also: case_large_ape
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

L  = 14.0;  % domain width (m)
H  = h_1+h_2;  % domain depth (m)
a_d = 0.019/2; % del_rho
delrho = abs(rho_1 - rho_2)/rho_2; %
a_d = delrho/2;
z0_d = h_1 + (h_pyc/2); % depth of pycnocline
d_d = h_pyc;
% Specify the general density profile which takes d_d as a second parameter
frho=@(z,d_d) 1-a_d*tanh((z+z0_d)/d_d);
frhoz=@(z,d_d) -(a_d/d_d)*sech((z+z0_d)/d_d).^2; % derivative in z

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%% Set up a range of sample APEs and calculate corresponding amplitudes
Ai = (1:10).*[3e-4; 0.5e-5; .5e-4]; % Set up an array of test APEs (m^4/s^2)
Ai = unique(Ai(:)); % Double check no duplicates & sort
verbose = 0;
calc_wavelength = false;
NX = 128; 
NZ = 64;
for i = 1:length(Ai) % Test for various APEs
    A = Ai(i);% APE for wave (m^4/s^2)
    
    %%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = clock;
    % Density profile for this wave with specified d_d
    rho  = @(z) frho(z, d_d);
    rhoz = @(z) frhoz(z, d_d);
        
    % Iterate the DJL solution
    djles_refine_solution;
    %djles_diagnostics; djles_plot; % uncomment to view progress at each step
       
    end_time = clock;
    fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));
    
    % Save loopy things
    amps(i, 1) = -wave_ampl;
end

%clearvars -except a_d Ai amps d_d delrho NX NZ fr* H h_* L m_path rho rhoz Ub* verbose wave_amp z0_d

%% Now run for this amplitude
ActualAPE = interp1(amps, Ai, wave_amp);
if isnan(ActualAPE)
    plot(amps, Ai, 'k-'); xlabel('Amplitude'); ylabel('APE');
    error('Amplitude not within tested range of APEs')
end
A = ActualAPE; % APE for wave (m^4/s^2)
calc_wavelength = true;
epsilon = 1e-4;
NX = 128; NZ = 64;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;

% Specify resolution and pycnocline parameter d_d according to this schedule
% Iterate the DJL solution
djles_refine_solution;
% Increase the resolution, reduce epsilon, and iterate to convergence
epsilon = 1e-5;
NX = 512; NZ = 128;
djles_refine_solution;

% Increase the resolution, and iterate to convergence
NX = 2048; NZ = 256;
djles_refine_solution;

% Compute and plot diagnostics
djles_diagnostics; djles_plot; 
wavelength = djles_wavelength(eta,L);

end_time = clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Save some loopy things

DJL.TestAPEs = Ai; 
DJL.TestAmps = amps;
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


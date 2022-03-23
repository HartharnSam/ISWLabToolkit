function [wavelength] = calc_DJL(h1, h2, h_pyc, amp, rho_1, rho_2)
%calc_DJL - One line description of what the function or script performs (H1 line)%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
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
% See also: case_sharp_pycnocline,  OTHER_FUNCTION_NAME2
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
addpath([m_path, '../08_DJL'])

%% Set up a range of sample APEs
Ai = (1:10).*[1e-4; 1e-5; .5e-4];
Ai = Ai(:);
lab_params; % Load in .m file parameters
verbose = 0;
for i = 1:length(Ai)
    A = Ai(i);% APE for wave (m^4/s^2)
    L  = 7.0;  % domain width (m)
    H  = h_1+h_2;  % domain depth (m)
    calc_wavelength = false;
    relax = 0.15; % use strong underrelaxation
    
    % Specify the general density profile which takes d_d as a second parameter
    a_d = 0.019/2; % del_rho
    a_d = delrho/2;
    z0_d = h_1; % depth of pycnocline
    
    frho=@(z,d_d) 1-a_d*tanh((z+z0_d)/d_d);
    frhoz=@(z,d_d) -(a_d/d_d)*sech((z+z0_d)/d_d).^2; % derivative in z
    
    % The velocity profile (zero for this case) (m/s)
    Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;
    
    %%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_time = clock;
    
    % Specify resolution and pycnocline parameter d_d according to this schedule
    NXlist = [  64   128    256     512];
    NZlist = [  32    64    128     256];
    ddlist = [0.01  0.01    h_pyc   h_pyc];
    for ddindex=1:length(ddlist)
        % Resolution for this wave
        NX = NXlist(ddindex);
        NZ = NZlist(ddindex);
        
        % Density profile for this wave with specified d_d
        d_d  = ddlist(ddindex);
        rho  = @(z) frho(z, d_d);
        rhoz = @(z) frhoz(z, d_d);
        
        % Iterate the DJL solution
        djles_refine_solution;
        %djles_diagnostics; djles_plot; % uncomment to view progress at each step
    end
    
    % Reduce epsilon, iterate to convergence
    epsilon=1e-5;
    djles_refine_solution;
    
    end_time=clock;
    fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));
    
    % Compute and plot the diagnostics
    %djles_diagnostics
    %djles_plot
    % Save some loopy things
    amps(i, 1) = -wave_ampl;
end

%% Now run for this amplitude
ActualAPE = interp1(amps, Ai, wave_amp);
A = ActualAPE;% APE for wave (m^4/s^2)
calc_wavelength = true;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;

% Specify resolution and pycnocline parameter d_d according to this schedule
NXlist=[  64    128     256     512];
NZlist=[  32    64      128     256];
ddlist=[0.01    0.01    h_pyc   h_pyc];
for ddindex=1:length(ddlist)
    % Resolution for this wave
    NX = NXlist(ddindex);
    NZ = NZlist(ddindex);
    
    % Density profile for this wave with specified d_d
    d_d  = ddlist(ddindex);
    rho  = @(z) frho(z, d_d);
    rhoz = @(z) frhoz(z, d_d);
    
    % Iterate the DJL solution
    djles_refine_solution
    %djles_diagnostics; djles_plot; % uncomment to view progress at each step
end

% Reduce epsilon, iterate to convergence
epsilon = 1e-5;
djles_refine_solution

% Raise resolution, iterate to convergence
NX = 2048; NZ = 256;
djles_refine_solution

end_time = clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics;
djles_plot;
% Save some loopy things
APE = A;
CalcAmp = -wave_ampl;
wavelength = djles_wavelength(eta,L);
max_u = max(u(:));

save('DJL', 'Ai', 'amps', 'APE', 'c', 'density','eta', 'etax', 'etaz', 'u', 'ux', 'uxze', 'uz', 'vorticity', 'w', 'wave_amp', 'wavelength', 'wx', 'wz','x', 'xc', 'z', 'zc');

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
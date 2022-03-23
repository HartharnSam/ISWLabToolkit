function [grd, e_dir, II, ok_comp] = calculate_dissipation(piv_dir, piv_fname, image_dir, image_fname, win_size, k_fit_range, nu, ice_props, diagnostic)
%CALCULATE_DISSIPATION - Calculates TKE (x, y) planes and plots images
%Using Doran et al., 2001 method to calculate firstly direct dissipation
%rate estimation, and then spectral estimate of the same. Carries out
%pre-processing of input PIV .dfi file, to remove edges, identify the
%pycnocline and ice edge.
%
% INPUTS:
%   piv_file - Single .dfi file
%   image_file - A corresponding raw image file .dfi
%   diagnostic - switch to plot, or not plot diagnostics
%   parameters - ordered list of the plots to make
%   x_lim, y_lim - x and y limits
%   im_ice_threshold - Threshold for brightness attributed to "Ice" (value
%                       out of 255)
%   ice_thickness - % Max thickness of ice from surface
%   win_size -  Size of dissipation calculation window (in pixels)
%   k_fit_range - Range [rad/m] over which to fit energy spectra for TKE dissipation calculation
%   nu - Kinematic viscosity
%
% Other m-files required: dfireadvel, dfi_grid_read, cmocean,
% find_boundaries, dissipation_gradient_2D, subaxis, spec2_ps_nopad,
% spec_ps_nopad, NasmythFit, subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Peter Sutherland, style revised Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 20-Sep-2021; Last revision: 20-Sep-2021
% MATLAB Version: 9.4.0.813654 (R2018a)
% -----------------------------
%% Set any defaults
% -----------------------------
if nargin < 5 
    win_size = 64; %Size of dissipation's PIV window
end
if nargin < 6
    k_fit_range = [300 800]; %Range [rad/m] over which to fit energy spectra for TKE dissipation alculation
end
if nargin < 7
    nu = 1.7E-6; % Kinematic viscosity
end

if nargin < 9 
    diagnostic = false;
end

%--------------------------------------------------
%%  END USER INPUTS
% -------------------------------------------------
%% Load in Velocity / Vorticity data
im = dfireadvel(fullfile(piv_dir, piv_fname)); %File name, read in uncompressed AND uncompacted piv file
U = im.cdata(:,:,1);
V = im.cdata(:,:,2);
%vort = im.cdata(:,:,3);

% Set world coordinates
[x_range, y_range, dx, dy, nx, ny, x0, y0, x, y]  = dfi_grid_read(im);
Grid = struct('x0', x0, 'y0', y0, 'dx', dx, 'dy', dy, 'nx', nx, 'ny', ny, ...
    'x_range', x_range, 'y_range', y_range);

% Re-allign matrix to correct orientation
if dx > 0
    %vort = flip(vort, 2);
    U = flip(U, 2);
    V = flip(V, 2);
end
if dy < 0
    %vort = flip(vort, 1);
    U = flip(U, 1);
    V = flip(V, 1);
end

% Blank off edge effects
U(1:3,:) = NaN; U((end-2):end,:) = NaN; U(:,1:3) = NaN; U(:,(end-2):end) = NaN;
V(1:3,:) = NaN; V((end-2):end,:) = NaN; V(:,1:3) = NaN; V(:,(end-2):end) = NaN;
%vort(1:3,:) = NaN; vort((end-2):end,:) = NaN; vort(:,1:3) = NaN; vort(:,(end-2):end) = NaN;
%clear vort

[m,n] = size(U);


%% Make the windowed grid
grd.x_pix = 3 + (win_size/2):(win_size/4):(n-win_size/2); % overlapping grid, offset from edge by 3
grd.y_pix = 3 + (win_size/2):(win_size/4):(m-win_size/2);
[grd.X_pix, grd.Y_pix] = meshgrid(grd.x_pix, grd.y_pix);
[grd.m, grd.n] = size(grd.X_pix);

grd.x = interp1([1, n], x_range, grd.x_pix);
grd.y = interp1([1, m], y_range, grd.y_pix);

% ----------------------------------------------------
%% Calculate Dissipation
% ----------------------------------------------------
%% ---- Direct dissipation estimate from gradients ---
e_dir.x = x; e_dir.y = y;
e_dir.e_filt = NaN(grd.m,grd.n);

[e_dir.e_diff, e_dir.e, e_dir.e_1d_x, e_dir.e_1d_y] = dissipation_gradient_2D(U, V, dx, dy, nu);

%% ---- Spectral dissipation estimates ----
grd.e_spec2D = NaN(grd.m, grd.n);
grd.e_spec2D_Nasmyth = NaN(grd.m, grd.n);

grd.eps_ux = grd.e_spec2D; grd.eps_uy = grd.e_spec2D;
grd.eps_vx = grd.e_spec2D; grd.eps_vy = grd.e_spec2D;
grd.eps_ux_Nasmyth = grd.eps_ux;
grd.eps_vy_Nasmyth = grd.eps_ux;
temp_dx = abs(dx); % TODO - check if this has any rationale...
temp_dy = abs(dy);
for lm = 1:grd.m
    for in = 1:grd.n
        M_temp = grd.y_pix(lm) + ((-win_size/2 +1):win_size/2);
        N_temp = grd.x_pix(in) + ((-win_size/2 +1):win_size/2);
        
        % 2D Spectra
        [S2u_t.S,S2u_t.kx,S2u_t.ky] = spec2_ps_nopad(U(M_temp,N_temp),2*pi/temp_dx,2*pi/temp_dy,win_size,win_size,'hanning');
        S2u_t.E11 = trapz(S2u_t.ky,S2u_t.S,1);
        S2u_t.k = S2u_t.kx(S2u_t.kx>0);
        S2u_t.E11 = 2*S2u_t.E11(S2u_t.kx>0);
        
        [S2v_t.S,S2v_t.kx,S2v_t.ky] = spec2_ps_nopad(V(M_temp,N_temp),2*pi/temp_dx,2*pi/temp_dy,win_size,win_size,'hanning');
        S2v_t.E11 = trapz(S2v_t.kx,S2v_t.S,2)';
        S2v_t.k = S2v_t.ky(S2v_t.ky>0);
        S2v_t.E11 = 2*S2v_t.E11(S2v_t.ky>0);
        
        % 1D spectra
        [Sx_t.k,Sx_t.Su,Sx_t.err_low_hi] = spec_ps_nopad(U(M_temp,N_temp)',2*pi/temp_dx,1,'hanning');
        Sx_t.Eu11 = mean(Sx_t.Su,2, 'omitnan');
        
        [~,Sx_t.Sv,~] = spec_ps_nopad(V(M_temp,N_temp)',2*pi/temp_dx,1,'hanning');
        Sx_t.Ev22 = mean(Sx_t.Sv,2, 'omitnan');
        
        [Sy_t.k,Sy_t.Su,Sy_t.err_low_hi] = spec_ps_nopad(U(M_temp,N_temp),2*pi/temp_dy,1,'hanning');
        Sy_t.Eu22 = mean(Sy_t.Su,2, 'omitnan');
        
        [~,Sy_t.Sv,~] = spec_ps_nopad(V(M_temp,N_temp),2*pi/temp_dy,1,'hanning');
        Sy_t.Ev11 = mean(Sy_t.Sv,2, 'omitnan');
        
        
        % Direct dissipation
        e_dir.e_filt(lm,in) = mean(mean(e_dir.e(M_temp,N_temp), 'omitnan'), 'omitnan');
        
        
        % k_fit_range = [200 900];
        ok_ku = (S2u_t.k>=min(k_fit_range))&(S2u_t.k<=max(k_fit_range));
        ok_kv = (S2v_t.k>=min(k_fit_range))&(S2v_t.k<=max(k_fit_range));
        
        % Fit k^{-5/3} curves to horizontal / vertical energy spectra
        grd.e_spec2D(lm,in) = real((mean((S2u_t.E11(ok_ku).*S2u_t.k(ok_ku).^(5/3)/(1.5*18/55)).^(3/2)) + ...
            mean((S2v_t.E11(ok_kv).*S2v_t.k(ok_kv).^(5/3)/(1.5*18/55)).^(3/2))) / 2 );
        
        % Fit Nasmyth spectrum to vertically / horizontally integrated 2D
        % spectra
        grd.e_spec2D_Nasmyth(lm,in) = (NasmythFit(S2u_t.k(ok_ku), S2u_t.E11(ok_ku), nu) + ...
            NasmythFit(S2v_t.k(ok_kv),S2v_t.E11(ok_kv), nu) ) / 2;
        % Calculates the best fitting epsilon for u and v spectra and
        % averages between them
        
        ok_kx = (Sx_t.k>=min(k_fit_range))&(Sx_t.k<=max(k_fit_range));
        ok_ky = (Sy_t.k>=min(k_fit_range))&(Sy_t.k<=max(k_fit_range));
        
        % Fit k^{-5/3} curves to 1-D spectra
        grd.eps_ux_c(lm,in) = mean( (Sx_t.Eu11(ok_kx).*Sx_t.k(ok_kx).^(5/3)/(1.5*18/55)).^(3/2) );
        grd.eps_ux(lm,in) = (Sx_t.k(ok_kx).^(-5/3)\Sx_t.Eu11(ok_kx)/(1.5*18/55)).^(3/2);
        
        % grd.eps_vy(lm,in) = mean((Sy_t.Ev11(ok_ky).*Sy_t.k(ok_ky).^(5/3)/(1.5*18/55)).^(3/2));
        grd.eps_vy(lm,in) = (Sy_t.k(ok_kx).^(-5/3)\Sy_t.Ev11(ok_kx)/(1.5*18/55)).^(3/2);
        
        % grd.eps_uy(lm,in) = mean((Sy_t.Eu22(ok_kx).*Sy_t.k(ok_kx).^(5/3)/(1.5*24/55)).^(3/2));
        grd.eps_uy(lm,in) = (Sy_t.k(ok_kx).^(-5/3)\Sy_t.Eu22(ok_kx)/(1.5*18/55)).^(3/2);
        
        % grd.eps_vx(lm,in) = mean((Sx_t.Ev22(ok_ky).*Sx_t.k(ok_ky).^(5/3)/(1.5*24/55)).^(3/2));
        grd.eps_vx(lm,in) = (Sx_t.k(ok_kx).^(-5/3)\Sx_t.Ev22(ok_kx)/(1.5*18/55)).^(3/2);
        
        % Fit Nasmyth spectra to horizontal / vertical 1D spectra
        grd.eps_ux_Nasmyth(lm,in) = NasmythFit(Sx_t.k(ok_kx),Sx_t.Eu11(ok_kx));
        grd.eps_vy_Nasmyth(lm,in) = NasmythFit(Sy_t.k(ok_kx),Sy_t.Ev11(ok_kx));
    end
end

grd.e_spec_1d = (grd.eps_ux + grd.eps_uy + grd.eps_vx + grd.eps_vy) / 4;
grd.e_spec_1d_Nasmyth = (grd.eps_ux_Nasmyth + grd.eps_vy_Nasmyth) /2;

ok_comp = (grd.X_pix>=min(k_fit_range))&(grd.X_pix<=max(k_fit_range)) ...
    &(grd.Y_pix>=min(k_fit_range))&(grd.Y_pix<=max(k_fit_range));

disp('==================================================================');
disp(['Window size = ' num2str(win_size)]);


function [e_dir, e_spec, ok_comp] = calc_dissipation(piv_dir, piv_fname, win_size, nu)
%CALC_DISSIPATION - Calculates TKE (x, y) planes and plots images
%Using Doran et al., 2001 method to calculate firstly direct dissipation
%rate estimation, and then spectral estimate of the same. Carries out
%pre-processing of input PIV .dfi file, to remove edges, identify the
%pycnocline and ice edge.
%
% NOTE: running the spectral part of the code is slow, so only run with 2
% arguments if using this. 
%
% INPUTS:
%   piv_dir - Directory for piv file
%   piv_fname - Single .dfi filename
%   win_size -  Size of dissipation calculation window (in pixels) [Default
%   = 64]
%   nu - Kinematic viscosity
%
% OUTPUTS:
%   e_spec - Structure of spectral - calculated dissipation, and relevant grids
%   e_dir - Structure of direct - calculated dissipation, and relevant grids
%   ok_comp - boolean ok calculation matrix
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
if nargin < 3
    win_size = 64; %Size of dissipation's PIV window
end
if nargin < 4
    nu = 1.7E-6; % Kinematic viscosity
end
if nargin < 5
    k_fit_range = [300 800]; %Range [rad/m] over which to fit energy spectra for TKE dissipation alculation
end
%--------------------------------------------------
%%  END USER INPUTS
% -------------------------------------------------
%% Load in Velocity / Vorticity data
if isa(piv_fname, 'struct')
    im = piv_fname;
else
    im = dfireadvel(fullfile(piv_dir, piv_fname)); %File name, read in uncompressed AND uncompacted piv file
end
U = im.cdata(:,:,1);
V = im.cdata(:,:,2);
%vort = im.cdata(:,:,3);

% Set world coordinates
[x_range, y_range, dx, dy, ~, ~, ~, ~, x, y]  = dfi_grid_read(im);

% Re-allign matrix to correct orientation
if dx > 0
    %vort = flip(vort, 2);
    U = flip(U, 2);
    V = flip(V, 2);
end
if dy > 0
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

%% Calculate relevant frequencies/wavenumbers
%dx = abs(im.xWorldPerPixel); dy = abs(im.yWorldPerPixel);
%f_nyq = pi*dx;
%Lx = abs(im.nx*im.xWorldPerPixel);
%k_fit_range = [2*pi/Lx 1/f_nyq];

%% Make the windowed grid
e_spec.x_pix = 3 + (win_size/2):(win_size/4):(n-win_size/2); % overlapping grid, offset from edge by 3
e_spec.y_pix = 3 + (win_size/2):(win_size/4):(m-win_size/2);
[e_spec.X_pix, e_spec.Y_pix] = meshgrid(e_spec.x_pix, e_spec.y_pix);
[e_spec.m, e_spec.n] = size(e_spec.X_pix);

e_spec.x = interp1([1, n], x_range, e_spec.x_pix);
e_spec.y = interp1([1, m], y_range, e_spec.y_pix);

% ----------------------------------------------------
%% Calculate Dissipation
% ----------------------------------------------------
%% ---- Direct dissipation estimate from gradients ---
e_dir.x = x; e_dir.y = y;
[e_dir.e_diff, e_dir.e_image, e_dir.e_1d_x, e_dir.e_1d_y] = dissipation_gradient_2D(U, V, dx, dy, nu);
e_dir.e_filt = NaN(e_spec.m,e_spec.n);

%% ---- Spectral dissipation estimates ----
e_spec.e2D = NaN(e_spec.m, e_spec.n);
e_spec.e2D_Nasmyth = NaN(e_spec.m, e_spec.n);

e_spec.eps_ux = e_spec.e2D; e_spec.eps_uy = e_spec.e2D;
e_spec.eps_vx = e_spec.e2D; e_spec.eps_vy = e_spec.e2D;
e_spec.eps_ux_Nasmyth = e_spec.eps_ux;
e_spec.eps_vy_Nasmyth = e_spec.eps_ux;
temp_dx = abs(dx); % TODO - check if this has any rationale...
temp_dy = abs(dy);
for lm = 1:e_spec.m
    for in = 1:e_spec.n
        M_temp = e_spec.y_pix(lm) + ((-win_size/2 +1):win_size/2);
        N_temp = e_spec.x_pix(in) + ((-win_size/2 +1):win_size/2);

        % Direct dissipation with window filter
        e_dir.e_filt(lm,in) = mean(mean(e_dir.e_image(M_temp,N_temp), 'omitnan'), 'omitnan');

        if nargout >= 2
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

            % k_fit_range = [200 900];
            ok_ku = (S2u_t.k>=min(k_fit_range))&(S2u_t.k<=max(k_fit_range));
            ok_kv = (S2v_t.k>=min(k_fit_range))&(S2v_t.k<=max(k_fit_range));

            % Fit k^{-5/3} curves to horizontal / vertical energy spectra
            e_spec.e2D(lm,in) = real((mean((S2u_t.E11(ok_ku).*S2u_t.k(ok_ku).^(5/3)/(1.5*18/55)).^(3/2)) + ...
                mean((S2v_t.E11(ok_kv).*S2v_t.k(ok_kv).^(5/3)/(1.5*18/55)).^(3/2))) / 2 );

            % Fit Nasmyth spectrum to vertically / horizontally integrated 2D
            % spectra
            e_spec.e2D_Nasmyth(lm,in) = (NasmythFit(S2u_t.k(ok_ku), S2u_t.E11(ok_ku), nu) + ...
                NasmythFit(S2v_t.k(ok_kv),S2v_t.E11(ok_kv), nu) ) / 2;
            % Calculates the best fitting epsilon for u and v spectra and
            % averages between them

            ok_kx = (Sx_t.k>=min(k_fit_range))&(Sx_t.k<=max(k_fit_range));
            ok_ky = (Sy_t.k>=min(k_fit_range))&(Sy_t.k<=max(k_fit_range));

            % Fit k^{-5/3} curves to 1-D spectra
            e_spec.eps_ux_c(lm,in) = mean( (Sx_t.Eu11(ok_kx).*Sx_t.k(ok_kx).^(5/3)/(1.5*18/55)).^(3/2) );
            e_spec.eps_ux(lm,in) = (Sx_t.k(ok_kx).^(-5/3)\Sx_t.Eu11(ok_kx)/(1.5*18/55)).^(3/2);

            % grd.eps_vy(lm,in) = mean((Sy_t.Ev11(ok_ky).*Sy_t.k(ok_ky).^(5/3)/(1.5*18/55)).^(3/2));
            e_spec.eps_vy(lm,in) = (Sy_t.k(ok_kx).^(-5/3)\Sy_t.Ev11(ok_kx)/(1.5*18/55)).^(3/2);

            % grd.eps_uy(lm,in) = mean((Sy_t.Eu22(ok_kx).*Sy_t.k(ok_kx).^(5/3)/(1.5*24/55)).^(3/2));
            e_spec.eps_uy(lm,in) = (Sy_t.k(ok_kx).^(-5/3)\Sy_t.Eu22(ok_kx)/(1.5*18/55)).^(3/2);

            % grd.eps_vx(lm,in) = mean((Sx_t.Ev22(ok_ky).*Sx_t.k(ok_ky).^(5/3)/(1.5*24/55)).^(3/2));
            e_spec.eps_vx(lm,in) = (Sx_t.k(ok_kx).^(-5/3)\Sx_t.Ev22(ok_kx)/(1.5*18/55)).^(3/2);

            % Fit Nasmyth spectra to horizontal / vertical 1D spectra
            e_spec.eps_ux_Nasmyth(lm,in) = NasmythFit(Sx_t.k(ok_kx),Sx_t.Eu11(ok_kx));
            e_spec.eps_vy_Nasmyth(lm,in) = NasmythFit(Sy_t.k(ok_kx),Sy_t.Ev11(ok_kx));
        end
    end
end

if nargout >= 2
    e_spec.e1d = (e_spec.eps_ux + e_spec.eps_uy + e_spec.eps_vx + e_spec.eps_vy) / 4;
    e_spec.e1d_Nasmyth = (e_spec.eps_ux_Nasmyth + e_spec.eps_vy_Nasmyth) /2;

    ok_comp = (e_spec.X_pix>=min(k_fit_range))&(e_spec.X_pix<=max(k_fit_range)) ...
    &(e_spec.Y_pix>=min(k_fit_range))&(e_spec.Y_pix<=max(k_fit_range));

    disp('==================================================================');
    disp(['Window size = ' num2str(win_size)]);
end


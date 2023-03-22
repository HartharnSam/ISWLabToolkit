function KdV = calc_kdv(Hs, h_2, d, delrho, amp_perc)
%CALC_KDV - % pseudospectral solution of the vertical structure
% of linear internal waves:
%  lambda(-D^2+k^2 I)phi = N^2(z)phi
% on 0<z<H with phi(0)=phi(H)=0
% This script varies the total depth for a fixed stratification if Hs is a
% vector, or calculates full physical solutions if Hs is a single number
%
% Inputs:
%    Hs - Depth of water column
%    h_2 - Height above bed of pycnocline centre in deepwater limit (H = Hs(1))
%    d - Pycnocline halfwidth
%    delrho - density difference (normalised)
%    amp_perc - Percentage of water column occupied by wave
% Outputs:
%    KdV - Output structure of KdV coefficients, and for physical solution
%    implementation, the physical solutions of the KdV equation
%
% Other m-files required: cmocean, djles_diffmatrix, dfles_gradient,
% subaxis, 
% Subfunctions: none
% MAT-files required: none
%
% See also: longwave_depthchange, djles_initial_guess, djles_diagnostics
% Author: Marek Stastna & Sam Hartharn-Evans
% Department of Applied Mathematics, University of Waterloo
% email address: mmstastna@uwaterloo.ca
% GitHub: https://github.com/HartharnSam
% 22-Mar-2022; Last revision: 22-Mar-2022
% MATLAB Version: 9.12.0.1884302 (R2022a)
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

close all;
format long, format compact

% set this to 1 if you want live plots
PLOT_NOW = 0;
save = false;
% This is what you need to get the differentiation matrix and the grid
Hdeep = Hs(1);
L = Hdeep*100;
g = 9.81; % Assume the lab is on Earth
delrho = delrho*0.5;
z0deep = h_2;
%d = 0.01; % d - pycnocline halfthickness

% Calculate Grids & Spatial things
N = 200; NX = 1000;
dx=L/NX;
xc = (.5:NX)*dx;  % Cell-centre grid
ks = (pi/L) * [0:(NX-1) -NX:-1];   % Wavenumber space
targetgrid='interior';

% here are the physical parameters and the scaling to the computational
% domain [-1, 1]

%Hs = linspace(Hdeep,0.2*Hdeep,mylen);

% Preallocate
r10_1 = NaN(length(Hs), 1); r01_1 = r10_1; r10_2 = r10_1; r01_2 = r10_1;
c1s = r10_1; c2s = r10_1; beta = r01_2; alpha = r01_2;

%% Run for each depth
for ii = 1:length(Hs)
    H = Hs(ii);
    % Re-do spatial domains in z for changing depths
    dz=H/N;
    zc = (0.5:N)*dz; [XC, ZC] = meshgrid(xc,zc);
    ms = (pi/H) * [0:(N-1) -N:-1]'; %
    zphys = zc(2:end-1)';

    % Differentiation matrices
    Dz   = djles_diffmatrix(dz, N, 1, 'not periodic');
    Dzz  = djles_diffmatrix(dz, N, 2, 'not periodic');
    Dzzc = Dzz(2:end-1, 2:end-1); % cut the endpoints for Dirichlet conditions

    % inline functions for the density and the derivative of the density
    z0 = z0deep+(H-Hdeep); % Pycnocline depth relative to current bed
    myrho=@(z) 1-delrho*tanh((z-z0)/d);
    myn2=@(z) (g*delrho/d)*sech((z-z0)/d).^2;

    if ii == 1
        figure
        plot(myrho(zphys), zphys);
    end

    n2physical=myn2(zphys);

    %% Solve KdV
    % make up the matrices for the e-val prog.
    %define B
    %define A
    B0 = sparse(diag(n2physical));
    B1 = sparse(B0.*0);
    B2 = sparse(Dzzc);
    % Solve the e-val prob
    [V, cc] = polyeig(B0, B1, B2);

    [c, ci]=sort(cc,'descend'); c_lw_1 = c(1); c_lw_2 = c(2);

    % Add boundary conditions
    phi1 = [0; V(:, ci(1)); 0];
    phi2 = [0; V(:, ci(2)); 0];

    % Normalise & compute E1
    E1 = abs(phi1)/max(abs(phi1));
    E2 = abs(phi2)/max(abs(phi2));

    if PLOT_NOW==1 % Plot the solution to the eigenvalue problem
        figure
        subplot(1,2,1)
        plot(n2physical,zphys+(Hdeep-H)),grid on,ylabel('z (m)'),xlabel('N^2')
        subplot(1,2,2)
        plot(phi1,zphys(2:end-1)+(Hdeep-H)),grid on,ylabel('z (m)'),xlabel('\phi')
        title(['c1 = ' num2str(c_lw_1, 4)])
    end

    % Here is WNL stuff calculated for each solution
    E1p = Dz*E1; E1p2 = E1p.^2; E1p3 = E1p.^3;
    E2p = Dz*E2; E2p2 = E2p.^2; E2p3 = E2p.^3;

    bot = sum((c_lw_1).*E1p2);
    r10_1(ii) = (-0.75/c_lw_1)*sum((c_lw_1).*(c_lw_1).*E1p3)/bot;
    r01_1(ii) = -0.5*sum((c_lw_1).*(c_lw_1).*E1.*E1)/bot;
    %fprintf('WNL gives: c_lw = %f, r10 = %f, r01 = %f\n\n', c_lw_1, r10_1, r01_1);

    bot = sum((c_lw_2).*E2p2);
    r10_2(ii) = (-0.75/c_lw_2)*sum((c_lw_2).*(c_lw_2).*E2p3)/bot;
    r01_2(ii) = -0.5*sum((c_lw_2).*(c_lw_2).*E2.*E2)/bot;
    %fprintf('WNL gives: c_lw = %2.2f m/s, r10 = %f, r01 = %f\n\n', c_lw_2, r10_2, r01_2);

    c1s(ii) = c_lw_1;
    c2s(ii) = c_lw_2;

    beta(ii) = r01_1(ii);
    alpha(ii) = r10_1(ii).*c_lw_1;
end

KdV = struct('Depths', Hs, 'c1s', c1s, 'c2s', c2s, 'r10_2', r10_2, 'r10_1', r10_1, 'r01_1', r01_1, 'r01_2', r01_2, 'c_0', c1s(1), 'alpha', alpha, 'beta', beta);
c_lw = c1s(1);

if length(Hs) == 1
    %% Calculate some physical parameters
    % Now optimise the b0, lambda parameters
    E = repmat(E1, 1, NX); E(:,1)=0; E(:,end)=0;
    if nargin < 5
        amp_perc = 0.05; % Start b0 as 5% of domain height
    end
    b0 = sign(r10_1(ii))*amp_perc*H;  
    lambda = sqrt( -6*r01_1 / (c_lw_1 * r10_1(ii) * b0) );
    c_nl = (1+(2/3)*r10_1*b0)*c_lw_1;

    %fprintf('b0 (Amp) = %2.2f m, lambda = %2.2f m, c_isw = %2.2f m/s \n', b0, lambda, c_nl);

    %% Calculate the corresponding fields (eta, density, u, w)
    eta = -b0*E.*sech((XC-L/2)/lambda).^2;
    eta(1,:)=0; eta(:,1)=0; eta(:,end)=0; eta(end,:)=0;

    density = myrho(ZC-eta); % density is the displaced density profile

    c   = sqrt(g*H/lambda); % DSS2011 Eq 24
    [etax, etaz] = djles_gradient(eta, ks, ms, 'odd', 'odd', targetgrid);
    u = c*etaz; % From djles_diagnostics
    w = -c*etax;

    KdV.U = u; KdV.W = w; KdV.density = density; KdV.X = XC; KdV.Z = ZC;
    KdV.b0 = b0; KdV.c_nl = c_nl; KdV.lambda = lambda;

    %% Plot those fields
    figure
    subaxis(3, 1, 1);
    pcolor(XC, ZC, density); cmocean('dense'); c = colorbar; ylabel(c, 'density');
    hold on
    contour(XC, ZC, density, [1 1], '-k');
    ylabel('z (m)');

    subaxis(3, 1, 2);
    pcolor(XC, ZC, u); caxis([-1 1].*.1); cmocean('balance'); c = colorbar;
    ylabel(c, 'u (m/s)'); ylabel('z (m)');

    subaxis(3, 1, 3);
    pcolor(XC, ZC, w); caxis([-1 1].*0.05); cmocean('balance'); c = colorbar;
    ylabel(c, 'w (m/s)'); ylabel('z (m)'); xlabel('x (m)');

end

if save
    save('KdV', 'KdV');
end
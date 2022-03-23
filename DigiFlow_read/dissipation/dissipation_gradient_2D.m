function [e_diff,e_image,e_1dx,e_1dy] = dissipation_gradient_2D(U,V,dx,dy, nu)
% dissipation_gradient_2D.m - function to compute dissipation using
% velocity gradients on a 2-D plane.
%
% The derivation for this expression can be found in:
%   Doron, P.; Bertuccioli, L.; Katz, J. & Osborn, T. R. Turbulence
%   Characteristics and Dissipation Estimates in the Coastal Ocean Bottom
%   Boundary Layer from PIV Data Journal of Physical Oceanography, 2001,
%   31, 2108-2134   
% 
% The major assumptions are:
%   1) The sampling resolution is finer than the Kolmogorv scale
%           eta = (nu^3/epsilon)^(1/4)
%      (E.g Tennekes and Lumley, Eqn 1.5.11)
%   2) Cross-stream gradients have the same magnitude as the measured
%   in-plane cross gradients.
% Doron et al. note that these assumptions are weaker than the isotropy
% assumptions of the inertial subrange fit and 1-D gradient method (e_1dx
% and e_1dy in this function).
%
%
% Inputs (all required):
%    U - U(y, x) horizontal velocity field (m/s) 
%    V - V(y, x) vertical velocity field (m/s)
%    dx - change in horizontal distance (m) per pixel
%    dy - change in vertical axis (m) per pixel
%
% Outputs:
%    e_diff - \epsilon_d (direct) estimate of the dissipation rate)
%    e_image - Image e(y, x) of \epsilon ("direct" estimate of dissipation before
%       spatial averaging)
%    e_1dx - Dissipation rate for isotrpoic hongeneous turbulence
%    e_1dy - 
%
% USAGE:
%   e_diff = dissipation_gradient_2D(U,V,dx,dy)
%
%   All units SI.
%
% Peter Sutherland - 2011/06/30
% Last revision: 22-Jun-2021 ; Sam Hartharn-Evans
if nargin<5
    nu = 1.7E-6; % [m^2/s]
end

% dissipation using gradients
du_dx = conv2(U,[1/2 0 -1/2],'same')./dx;
du_dy = conv2(U,[1/2 0 -1/2]','same')./dy;
dv_dx = conv2(V,[1/2 0 -1/2],'same')./dx;
dv_dy = conv2(V,[1/2 0 -1/2]','same')./dy;

du_dx(:,1) = NaN; du_dx(:,end) = NaN; 
dv_dx(:,1) = NaN; dv_dx(:,end) = NaN; 
du_dy(1,:) = NaN; du_dy(end,:) = NaN; 
dv_dy(1,:) = NaN; dv_dy(end,:) = NaN; 

%e_diff = 3*nu*(mean(du_dx(:).^2, 'omitnan') + mean(du_dy(:).^2, 'omitnan') + mean(dv_dx(:).^2, 'omitnan') +...
%    mean(dv_dy(:).^2, 'omitnan') + mean(2*du_dy(:).*dv_dx(:), 'omitnan') + 2/3*mean(du_dx(:).*dv_dy(:), 'omitnan'));

% Version from Doran et al., 2001 - 2D by assuming isostophy
e_image = 3*nu*(du_dx.^2 + du_dy.^2 + dv_dx.^2 + dv_dy.^2 + 2*du_dy.*dv_dx + 2/3*du_dx.*dv_dy);
% Version from SPINS - 2D by removing any v or y terms
%e_image = 2*nu*(du_dx.^2 + dv_dy.^2 + 2*(.5*(du_dy + dv_dx)).^2);

e_diff = mean(e_image(:) , 'omitnan');

e_1dx = 15*nu*(du_dx.^2); % dissipation rate for homogeneous isotropic turbulence; 
e_1dy = 15*nu*(dv_dy.^2);


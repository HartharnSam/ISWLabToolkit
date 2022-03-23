function [nasmyth, universal, universal_orig] = NasmythSpectrum()
% NasmythSpectrum - A simple function to produce Nasmyth's universal spectrum.
% Originally shown graphically in Nasmyth 1970 and then in table
% form in Oakey 1982 (Table A1).
%
% Also creates a structure called universal which is the Nasmyth spectrum
% in the coordinates of Pope 2000.
%
% Outputs:
%   nasmyth - a structure which contains the data from Oakey's
% table.
%   universal - a structure which is the Nasmyth spectrum in the coordinates of Pope 2000.
%   universal_orig - as universal, but in "modern units"
%
% k_ks = \hat{k}/k_s
%
% ks = (epsilon*nu^(-3))^(1/4);
%
% % velocity spectrum
% phi_11 = (epsilon*nu^5)^(1/4)*F;
% phi_22 = (epsilon*nu^5)^(1/4)*F2;
%
% shear or gradient spectrum
% (2*pi*k)^2*phi_22 = ks^2*(epsilon*nu^5)^(1/4)*G2
%

nasmyth.kh_knh = [2.83E-4;
    5.033E-4;
    8.950E-4;
    1.592E-3;
    2.830E-3;
    5.033E-3;
    8.950E-3;
    1.592E-2;
    2.830E-2;
    5.033E-2;
    7.977E-2;
    1.264E-1;
    1.592E-1;
    2.004E-1;
    2.522E-1];


nasmyth.F = [1.254E+5;
    4.799E+4;
    1.842E+4;
    7.050E+3;
    2.699E+3;
    1.036E+3;
    3.964E+2;
    1.490E+2;
    3.574E+1;
    5.600E+0;
    7.214E-1;
    6.580E-2;
    1.812E-2;
    4.552E-3;
    1.197E-3];

nasmyth.F2 = [1.671E+5;
    6.397E+4;
    2.455E+4;
    9.404E+3;
    3.598E+3;
    1.380E+3;
    5.320E+2;
    2.302E+2;
    6.880E+1;
    1.373E+1;
    2.101E+0;
    2.127E-1;
    6.161E-2;
    1.570E-2;
    4.011E-3];

nasmyth.G2 = [0.5285;
    0.6397;
    0.7763;
    0.9404;
    1.138;
    1.380;
    1.682;
    2.302;
    2.176;
    1.373;
    0.5278;
    0.1342;
    0.0616;
    0.0249;
    0.0101];




% universal spectrum in "modern" units
universal_orig.kn = 2*pi*nasmyth.kh_knh;
universal_orig.F = nasmyth.F/2/pi;


% extend the universal spectrum inertial subrange
universal.kn = [logspace(-5,-3,10)'; universal_orig.kn];
universal.F = [universal_orig.F(1)*universal_orig.kn(1)^(5/3)*universal.kn(1:10).^(-5/3); universal_orig.F];



% % % Back to real units:
% doplot = false;
% if doplot
%
%     nu = 1E-6;                      % Viscosity
%     epsilon = 10.^(-8:-2);          % dissipation [W/kg]
%     eta = (nu^3./epsilon).^(1/4);   % Kolmogorov scale
%     k = universal.kn*(1./eta);      % wavenumber [rad/m]
%     E_11 = universal.F*(epsilon*nu^5).^(1/4);   % Energy spectrum
%     figure
%     loglog(k, E_11);
%     xlabel('Wavenumber (rad/m)')
%     ylabel('Energy')
%     for i = 1:length(epsilon)
%         leg_item{i} = ['$\epsilon = $', num2str(epsilon(i))];
%     end
%
%     legend(leg_item, 'interpreter', 'latex')
% end








function e = NasmythFit(k,E11, nu, isPlot)
% NASMYTHFIT - Calculates value of e that closest fits the the data
% to the NasmythSpectrum
% Outputs:
%    e - Dissipation
%
% Other m-files required: NasmythSpectrum
% Subfunctions: nd_spec_diff
% MAT-files required: none
%
% Author: Peter Sutherland

% Call Nasmyth Spectrum structure
if nargin < 4
    isPlot = false;
end
if nargin < 3
    nu = 1.7E-6; % Kinematic viscosity
end
[~, universal] = NasmythSpectrum;

% Stop iterating at 1000 iterations
opts = optimset('MaxIter',1000); % optimset(''Display','iter',TolFun',1E-10 )
e_lim = [-8 2];
[log10_e,~,exitflag] = fminbnd(@nd_spec_diff, e_lim(1), e_lim(2), opts, k, E11, universal, nu);

e = 10^log10_e;
if exitflag<0, error('Failed to converge'); end


if isPlot
    figure(99); clf;
    loglog(universal.kn,universal.F,'k');
    hold on;
    e_test = linspace(-10,1,100);
    fc = NaN(1:length(e_test));
    for i = 1:length(e_test)
        fc(i) = nd_spec_diff(e_test(i),k,E11,universal);
    end
    plot(e_test,(fc),'.-');
    
    nu = 1E-6;
    k_nd = k*(e*nu^(-3)).^(-1/4);
    E_nd = E11*(e*nu^5)^(-1/4);
    plot(k_nd,E_nd,'r--','LineWidth',2);
end

end


function sq_cost = nd_spec_diff(log10_e,k,E11,universal, nu)
% Calculates a difference between the spectral estimate and the observed
% for given e
%nu = 1E-6;

e = 10^log10_e;

k_nd = k*(e*nu^(-3)).^(-1/4);
E_nd = E11*(e*nu^5)^(-1/4);
% loglog(k_nd,E_nd);

F_tmp = exp(interp1(log(universal.kn),log(universal.F),log(k_nd), 'linear', 'extrap'));
% hold on
% loglog(k_nd,F_tmp);
sq_cost = log10(sum((F_tmp - E_nd).^2));
%pause(0.05);

end












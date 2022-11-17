% Calculate N2 from lab probe measurements
clc; clearvars; close all;

% CHANGE HERE the probe measurement ID
load 11102022

% Loads in the data from the .mat file
z = data.FittedData(:, 4);
rho = data.FittedData(:, 2);

% Initial plot of density
subplot(1, 2, 1)
plot(rho, z, '.k');

% Sort & smooth (aggressively)
%[z, zind] = unique(z, 'sorted');
%rho = rho(zind);
rho = smooth(z, rho, 0.10, 'rloess');

% Plot smoothed density profile
hold on
plot(rho, z, '-k');

% Calculate derivatives using the FiniteDiff tool (in SPINSmatlab)
Dmat = FiniteDiff(z, 1, 2, false, false);
drho_dz = Dmat*rho;
N2 = (-9.81/1026)*drho_dz;

% Plot the N2
subplot(1, 2, 2);
plot((N2), z);

fprintf('Max N^2 = %2.2f s^-2 \n', max(N2))

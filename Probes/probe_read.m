function [raw_data, fitted_data, outputs] = probe_read(input_filename, travel, total_depth, density_lims, cut, probes, depth_shift, pot_top_height, calibrated)
%PROBE_READ - Routine for comparing the initial density profiles of each probe.
%Function which takes in arguments from physical readings, to output the
%stratification in the tank. Uses hydrometer readings to calibrate the probes
%
% Syntax:  [raw_data, fitted_data, gridded_fit] = probe_read(input_filename, travel, total_depth, density_lims, cut, probes, depth_shift, pot_top_height, calibrated)
%
% Inputs:
%    input_filename - Name of file within current directory containing raw
%    probe data (in text/csv format)
%
%    travel - Total travel of potentiometer
%    total_depth - Total water depth
%    density_lims - [Lower layer density; upper layer density]
%    cut - Number of points to cut from [top bottom] of profile
%    probes - list of channels that probe data is on (e.g. [2 3])
%    depth_shift - some pointwise depth correction
%    pot_top_height - Probe height at top reading (usually = total_depth)
%    calibrated - [optional, default false] switch for using a known calibration of potentiometer
%    readings to height above bed
%
% Outputs:
%    raw_data - copy of unprocessed probe data
%    fitted_data - copy of the data in height above bed, density format
%    outputs - a structure containing information about the profile (e.g.
%    pycnocline thickness)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%---------------------------------------------------
%% BEGIN CODE %%
%--------------------------------------------------
% Load in probe data
inputs = importdata(input_filename);
if contains(input_filename, '.csv')
    inputs = inputs.data;
end
% Parse function inputs
if nargin < 9
    calibrated = false;
end

density_top = density_lims(1);
density_bot = density_lims(2);

%% Average top and bottom readings
% Number of points and channels in probe reading
Input = inputs;
%Input(:, 1) = sgolayfilt(Input(:, 1),1,5);

% Number of points and channels in probe reading
[N_points, ~] = size(Input); % For cut dataset
% later for top/bottom calibration)

% average the input parameters used for calibration,
% using multiple measurements from the top and bottom
N_avg = 100;
% average the top
Input_top_avg = mean(Input(1:N_avg, :));
%pot_top_avg =  mean(inputs([1:n_avg]+cut(1)-(n_avg/2), 1)); %(1);    % pot-meter reading at top
pot_top_avg =  Input_top_avg(1);    % pot-meter reading at top

% average the bottom

Input_bot_avg = mean(Input(N_points-N_avg+1:N_points, :));

% find the average difference between top and bottom readings
Input_avg_diff = Input_bot_avg - Input_top_avg;
pot_avg_diff = Input_avg_diff(1);   % pot-meter average difference

% replicate onto a matrix for the conversion
Input_top_avg_mat  = repmat(Input_top_avg, N_points,1);
Input_avg_diff_mat = repmat(Input_avg_diff,N_points,1);

%% convert readings to physical values based on hydrometers
% values are forced to be correct at the ends
den_diff = density_bot - density_top;
den_avg = density_top + den_diff/2;

% do conversion for all columns, then fix the pot-meter reading
phys      = density_top + (Input - Input_top_avg_mat)./Input_avg_diff_mat*den_diff;
if ~calibrated
    phys(:,1) = pot_top_height - (Input(:,1) - pot_top_avg)/pot_avg_diff*travel - depth_shift;
    phys(:, 1) = smooth(phys(:, 1), 10);
else
    % Old calibration - phys(:,1) = -((.3709.*(Input(:,1).^2))-(1.5815.*Input(:,1))+1.3831);
    %phys(:, 1) = ....
end

% remove points from top or bottom
cut_top = cut(1);
cut_bot = cut(2);
if cut_top > 0 || cut_bot > 0

    N_new = N_points - cut(1) - cut(2);       % New number of points
    new_top_avg = mean(phys(cut_top+1:cut_top+N_avg+1, :));                 % average of correct top section
    new_bot_avg = mean(phys(N_points-cut_bot-N_avg+1:N_points-cut_bot, :)); % average of correct bottom section
    new_avg_diff = new_bot_avg - new_top_avg;           % difference of the averages
    % replicate onto a matrix for the new conversion
    new_top_avg_mat  = repmat(new_top_avg, N_new,1);
    new_avg_diff_mat = repmat(new_avg_diff,N_new,1);

    % do new conversion to fix incorrect edges
    new_input = phys(cut_top+1:N_points-cut_bot, :);
    phys         = density_top + (new_input - new_top_avg_mat)./new_avg_diff_mat*den_diff;
    phys(:,1) = new_input(:,1); % return the correct depth readings into the matrix
    N_points = N_new;           % update number of points
end

%% Plot raw readings
figure(1)
clf, hold on
clr = ['b','k','r','g','y','c','m'];

ii=1;
for probe = probes
    plot(phys(:,probe), phys(:,1),'.','Color',clr(ii))
    ii=ii+1;
end
grid off,

set(gca, 'XDir', 'normal');
xlabel('Density (kg/m^3)'),
ylabel('HAB (m)'),
ylim([0 total_depth])
title('Initial stratification in the main tank')
leg = cellstr(num2str(probes'-1, 'probe %-d'));
legend(leg, 'location', 'southwest')

%% Fit to data
% set fitting region
fit_perc = 0.25;    % percentage of edge of profile to ignore in fit
top_edge = density_top + fit_perc*den_diff; % top edge of density fit
bot_edge = density_bot - fit_perc*den_diff; % bottom edge of density fit

% initialize vectors
ii=1;
fit_top = zeros(1,length(probes));
fit_bot = fit_top;
fit_mid = fit_top;
rho_mid = fit_top;

tanh_fit = zeros(2,length(probes));
% loop over all probes
for probe = probes
    % indices for top edge
    [~, top_ind] = min(abs(phys(:,probe) - top_edge));
    % indices for bottom pycnocline
    [~, bot_ind] = min(abs(phys(:,probe) - bot_edge));
    % list of indices to use in fit
    inds = top_ind:bot_ind;

    % do linear fits
    [p,~] = polyfit(phys(inds,probe), phys(inds,1), 1);
    % add to plot
    rho_line = linspace(density_top, den_avg, N_points);
    %plot(rho_line, p(1)*rho_line + p(2),'--', 'Color', clr(ii));

    % locations of edges of pycnoclines
    fit_top(ii) = p(1)*density_top + p(2);
    fit_bot(ii) = p(1)*den_avg + p(2);
    fit_mid(ii) = (fit_top(ii)+fit_bot(ii))/2;
    % location of pycnocline centres from density readings
    [~, rho_mid_ind] = min(abs(phys(:,probe) - den_avg));
    rho_mid(ii) = phys(rho_mid_ind,1);

    % do hyperbolic tangent fit (or fit with the error function)
    func = @(param, z) den_avg - den_diff/2*tanh((z-param(1))/param(2));
    %func = @(param, z) den_avg - den_diff/2*erf((z-param(1))/param(2));
    param = [rho_mid(ii) (rho_mid(ii)-fit_bot(ii))];    % initial guess
    tanh_lb = [0 0]; tanh_ub = [total_depth total_depth]; % Set physical bounds for the pycnocline width and depths as total water depth
    tanh_fit(:,ii) = lsqcurvefit(func, param, phys(:,1), phys(:,probe), tanh_lb, tanh_ub);
    plot(func(tanh_fit(:, ii), phys(:,4)), phys(:,4), ':', 'Color', clr(ii), ...
        'DisplayName',num2str(probes(ii)'-1, 'tanh %-d'));
    ii=ii+1;
end

%% Return useful information
Lower_Depth = p(1)*density_bot + p(2); %returns depth at which best fit line = bottom density
Upper_Depth = total_depth - (p(1)*density_top + p(2)); %returns depth at which best fit line = bottom density
Pycnocline_Thickness = total_depth - (Lower_Depth+ Upper_Depth);

% print out info from Linear Fit
fprintf('----- ---------- ------ \n');
fprintf('\n')
fprintf('Average density: %i kg/m^3 \n',den_avg)
fprintf('----- ---------- ------ \n');
fprintf('\n----- Linear Fit ------ \n');
fprintf('----- ---------- ------ \n');
fprintf('Lower Layer Thickness = %2.4f m \n', Lower_Depth);
fprintf('Upper Layer Thickness = %2.4f m \n', Upper_Depth);
fprintf('Pycnocline Thickness = %2.4f m \n', Pycnocline_Thickness);
% print out info from Tanh Fit
fprintf('----- ---------- ------ \n');
fprintf('\n----- Tanh Fit ------ \n');
fprintf('----- ---------- ------ \n');
fprintf('Pycnocline Depth = %2.4f m \n', tanh_fit(1));
fprintf('Pycnocline Thickness = %2.4f m \n', tanh_fit(2));
fprintf('----- ---------- ------ \n');

%% Function outputs

raw_data = Input;
fitted_data = phys;
outputs.LinearFit.LowerThick = Lower_Depth;
outputs.LinearFit.UpperThick = Lower_Depth;
outputs.LinearFit.PycThick = Pycnocline_Thickness;
outputs.TanhFit.pyc_depth = tanh_fit(1);
outputs.TanhFit.pyc_thick = tanh_fit(2);

%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
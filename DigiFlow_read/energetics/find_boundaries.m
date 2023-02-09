function II = find_boundaries(image_fname, type, ice_thickness, im_ice_threshold, diagnostic)
%FIND_BOUNDARIES - Using brightness of the raw image, identifies the location
% of the pycnocline or ice edge
%
% Syntax:  II = find_boundaries(image_fname, type, ice_thickness, im_ice_threshold, diagnostic)
%
% Inputs:
%    image_fname - Filename of raw image (String)
%    type - "ice" or "pycnocline" - boundary to be detected
%    ice_thickness - greatest depth (from top) where ice could be antipated
%    [default = 0.5]
%    im_ice_threshold - brightness threshold value for ice [default = 100]
%    diagnostic -  switch for image diagnostics
%
% Outputs:
%    II - Structure containing:
%           - II.og - dfireadvel structure
%           - II.im - dfireadvel colordata (im.cdata)
%           - II.im - dfireadvel colordata (im.cdata)
%           - II. grid - outputs from dfi_grid_read
%           - II.y_ice or II.y_interface - Location of detected interface
%           - II.qpy / II.qpx - Region bounded by ice/upper boundary
%           - II.y_shallow_peak or II.y_deep_peak - image used for peak
%           detection
%
% Other m-files required: dfireadvel, cmocean, newbluewhitered, plot_dfi,
% subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Peter Sutherland; Re-styled Sam Hartharn-Evans
%
%
% 21-Sep-2021; Last revision: 21-Sep-2021
% MATLAB Version: 9.4.0.813654 (R2018a)
%
%% Check args in

if nargin < 3 || (isempty(ice_thickness) && strcmpi(type, 'ice'))
    ice_thickness = .05;
end
if nargin < 4 || (isempty(im_ice_threshold) && strcmpi(type, 'ice'))
    ice_thickness = 100;
end
if nargin < 5
    diagnostic = false;
end

%% Read in the "image"
if ~isa(image_fname, 'struct')
    II.og = dfireadvel(fullfile(image_fname));
    II.im = II.og.cdata;
    [II.x_range, II.y_range, II.dx, II.dy, II.nx, II.ny, II.x0, II.y0, II.x, II.y] = dfi_grid_read(II.og);

else
    II = image_fname;
end
[II.m,II.n] = size(II.im);

% Re-allign matrix to correct orientation (important orientation is correct
% for identifying interfaces etc. later on)
if II.dx < 0
    II.im = flip(II.im, 2);
end
if II.dy < 0
    II.im = flip(II.im, 1);
end

%% Find interfaces using image brightness
locs_s = NaN(II.nx, 1); loc_thr = nan(II.nx, 1); locs_d = NaN(II.nx, 1);

tmp.image_peak = II.im; % Will be the ice
tmp.y_interface = NaN(II.nx, 1);

if license('test', 'image_toolbox')
    tmp.image_peak = imgaussfilt(II.im, 5, 'FilterSize', 11); % filter to smooth image (and thus the detected edges)
end

switch lower(type)
    case 'ice'
        
        ice_depth = max(II.y_range)- ice_thickness;
        tmp.image_peak(II.y<ice_depth, :) = NaN; % Remove data away from surface
        for i_pix = 1:II.nx
            % ice (brightness threshold from top of tank)
            try
                loc_thr(i_pix) = find(tmp.image_peak(:, i_pix) > im_ice_threshold, 1, 'last');
                tmp.y_interface(i_pix) = II.y(loc_thr(i_pix));
            catch
            end
        end
        II.qpy = [tmp.y_interface ones(II.nx, 1)*max(II.y_range)]'; % Identifies region filled by ice
        II.image_shallow_peak = tmp.image_peak;
        II.y_ice = tmp.y_interface;
        
    case 'pycnocline' % pycnocline (brightness peak)
        tmp.image_peak(II.y > max(II.y)*.9 | II.y < max(II.y)*.1, :) = NaN; % Remove data near-surface + bed
        for i_pix = 1:II.nx
            try
                [~, locs_d(i_pix)] = findpeaks(double(tmp.image_peak(:, i_pix)), 'SortStr', 'descend', 'NPeaks', 1);
                tmp.y_interface(i_pix) = II.y(locs_d(i_pix));% Y Location of the lower brightness peak (pycnocline)
            catch
            end
        end
        II.image_deep_peak = tmp.image_peak;
        II.y_interface = tmp.y_interface;
        
end

%% Select ice bottom and layer interface based on ice thickness

II.qpx = [1; 1]*II.x; % Paired X coordinates at each x point

%% Plot Raw Image with boundaries marked
if diagnostic
    figure(11);
    % First plot the full raw image with boundaries on
    diag_sp(3) = subaxis(2, 2, 3);
    imagesc(II.x, II.y, II.im); shading flat; cmocean('grey');
    set(gca, 'YDir', 'normal', 'XDir', 'reverse');
    %axis equal;
    xlim([min(II.x_range) max(II.x_range)]); ylim([min(II.y_range) max(II.y_range)]);
    hold on
    
    switch lower(type)
        case 'ice'
            plot(II.x, II.y_ice,'w','LineWidth',4,'HandleVisibility','Off');
            p(2) = plot(II.x,II.y_ice,'b','LineWidth',1);
            legend([p(2)], 'Ice Bottom', 'Location','w')
            xlabel('x [m]');
            ylabel('z [m]');
            title('Detected peaks')
            
            % Plot the shallow peak fitting image
            diag_sp(1) = subaxis(2, 2, 1);
            imagesc(II.x, II.y, II.image_shallow_peak); shading flat; cmocean('grey');
            set(gca, 'YDir', 'normal', 'XDir', 'reverse');
            axis equal;
            xlim([min(II.x_range) max(II.x_range)]); ylim([min(II.y_range) max(II.y_range)]);
            xlabel('x [m]');
            ylabel('z [m]');
            title('Shallow Peak (ice detection) Image')
            
            
        case 'pycnocline'
            plot(II.x,II.y_interface,'w-','LineWidth',3,'HandleVisibility','Off');
            p(1) = plot(II.x,II.y_interface,'k--','LineWidth',1);
            legend(p(1),'Pycnocline','Location','w')
            xlabel('x [m]');
            ylabel('z [m]');
            title('Detected peaks')
            % Plot the deep peak fitting image
            diag_sp(2) = subaxis(2, 2, 2);
            imagesc(II.x, II.y, II.image_deep_peak); shading flat; cmocean('grey');
            set(gca, 'YDir', 'normal', 'XDir', 'reverse');
            %axis equal;
            xlim([min(II.x_range) max(II.x_range)]); ylim([min(II.y_range) max(II.y_range)]);
            xlabel('x [m]');
            ylabel('z [m]');
            title('Deep Peak (pyc. detection) Image')
    end
    
end
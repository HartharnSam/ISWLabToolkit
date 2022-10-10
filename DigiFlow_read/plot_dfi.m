function [ax, BG_data, grid, Stats] = plot_dfi(fnm_in, BackgroundType, ForegroundType, ArgsIn, grid,ax )
%PLOT_DFI - Reads in and plots DigiFlow Floating Point Image (.dfi) Data. Gives plot
% on defined axes, can handle PIV, timeseries or image files
%
% Inputs:
%    fnm_in - Filename of data, can be {'Background data', 'Foreground data'}
%    BackgroundType - Creating Image Process for background ('TimeSeries'
%    or 'Image' or 'Vorty', or 'U', or 'V', or 'Speed', or 'Dissipation' or
%    'none')
%    ForegroundType - Foreground process: 'quiver', 'streamline', 'none'  [REQUIRED]
%    ArgsIn - Structure containing additional inputs required for each file
%    grid - [OPTIONAL] - Use pre-existing form of the "Grid" (for speed)
%    ax - [OPTIONAL] - use axes handle
%
% Outputs:
%    ax - Axes handle data plotted on
%	 data - Raw data relating to the "parameter" (for image and timeseries, this is the only data)
%	 grid - Structure containing spatial reference information for the data
%    Stats - Structure containing some statistics on the output data
%
%
% Other m-files required: dfireadvel, dfi_grid_read
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 05-Nov-2020; Last revision: 20-Jan-2022
% MATLAB Version: 9.9.0.1467703 (R2020b)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%% Initialise Script
if strcmpi(ForegroundType, "quiver") || strcmpi(ForegroundType, "streamline")
    q_scale = ArgsIn.q_scale; % Size of quivers
    dn = ArgsIn.Resolution; % quiver resolution
else
    dn = 1;
end

if nargin < 6
    ax = gca;
end
axes(ax)
%% Read in File
im = cell(length(fnm_in), 1);
%x = cell(length(fnm_in), 1); y = cell(length(fnm_in), 1);

for ii = 1:length(fnm_in)
    if iscell(fnm_in)
        im{ii} = dfireadvel(fnm_in{ii});  %File name, read in uncompressed AND uncompacted dfi file
    else
        im{ii} = dfireadvel(fnm_in);
    end
    
    %% Set World Coordinate System (if not input as variable)
    if ~exist('grid', 'var') || length(grid) < ii || isempty(grid{ii})
        try
            grid{ii}  = dfi_grid_read(im{ii});
            X6 = grid{ii}.X; Y6 = grid{ii}.Y;
            x26 = X6(1:dn:size(X6,1), 1:dn:size(X6,2));
            y26 = Y6(1:dn:size(X6,1), 1:dn:size(X6,2));
            grid{ii}.x2 = x26; grid{ii}.y2 = y26; clear X6 Y6 x26 y26
            x{ii} = grid{ii}.x; y{ii} = grid{ii}.y;
            X{ii} = grid{ii}.X; Y{ii} = grid{ii}.Y;
            x2{ii} = grid{ii}.x2; y2{ii} = grid{ii}.y2;
            
        catch
            warning('No grid to read in image data')
        end
        
        %% Begin Script
        if (ii == 1 && ~strcmp(BackgroundType, 'none'))% Check background type
            switch lower(BackgroundType)
                % 'TimeSeries' or 'Image' or 'Vorty', or 'U', or 'V', or 'Speed', or 'Dissipation' or 'none'
                case "timeseries"
                    try
                        parameter = ArgsIn.Parameter; % Background plotting parameter
                    catch
                        parameter = 1;
                    end
                    [~, search_fnm] = fileparts(fnm_in{ii});
                    charstring = char(strcat('Output[\d] := "', search_fnm));
                    output_no = im{ii}.CreatingProcess(regexpi(im{ii}.CreatingProcess, charstring)+6);
                    
                    % Identify time direction
                    charstring = char(strcat('TimeDirection', convertCharsToStrings(output_no), ' := "'));
                    [~, timeDirection_ind] = regexpi(im{ii}.CreatingProcess, charstring);
                    time_direction = im{ii}.CreatingProcess(timeDirection_ind+(1:4));
                    
                    % Identify timeseries type (row, column, y(s) = F(x))
                    charstring  = char(strcat('Kind', convertCharsToStrings(output_no), ' := "'));
                    [~, kind_ind] = regexpi(im{ii}.CreatingProcess, charstring);
                    kind = im{ii}.CreatingProcess(kind_ind+(1:3));
                    switch lower(kind)
                        case 'y(x'
                            disp('Slope Timeseries')
                            %% Need to reset WCS
                            [~, xMaxInd] = regexpi(im{ii}.CreatingProcess, ...
                                char(strcat('xMax', convertCharsToStrings(output_no), ' :=')));
                            xMax = str2double(im{ii}.CreatingProcess(xMaxInd+(1:14)));
                            [~, xMinInd] = regexpi(im{ii}.CreatingProcess, ...
                                char(strcat('xMin', convertCharsToStrings(output_no), ' :=')));
                            xMin = str2double(im{ii}.CreatingProcess(xMinInd+(1:14)));
                            
                            % Collects the y=mx+c equation from the initial dialog response
                            [y_of_x_1, y_of_x_2] = regexpi(im{ii}.CreatingProcess, ...
                                char(strcat('y_of_x', convertCharsToStrings(output_no), ' := "(x*[()0123456789.y])')));
                            %                 %char(strcat('y_of_x', convertCharsToStrings(output_no), ' := "[()0123456789.y]+*x')) %|
                            inds = (y_of_x_1+15):(y_of_x_2+9);
                            %i = regexpi(im.CreatingProcess((y_of_x_1+12):(y_of_x_2-2)), ['\d']);
                            i = regexpi(im{ii}.CreatingProcess(inds), '\d');
                            
                            inds = inds(1:max(i));
                            x_gradient = str2double(im{ii}.CreatingProcess(inds));%[min(i):max(i)]));
                            
                            dist = [0 sqrt((xMax-xMin).^2+((xMax-xMin)*x_gradient).^2)];
                            switch time_direction
                                case "To r"
                                    y = dist;
                                    flip(x);
                                case "Upwa"
                                    x = dist;
                                otherwise
                                    error("no time direction set")
                            end
                        case {"row", "col"}
                            disp('Column or Row Timeseries')
                    end
                    
                case "image"
                    disp('Image File')
                    try
                        parameter(ii) = ArgsIn.Parameter; % Background plotting parameter
                    catch
                        parameter(ii) = 1;
                    end
                    grid{ii}.FrameRate = 1/im{ii}.tStep;
                case "u"
                    parameter(ii) = 1;
                    grid{ii}.FrameRate = 1/im{ii}.tStep;
                    
                case "v"
                    parameter(ii) = 2;
                    grid{ii}.FrameRate = 1/im{ii}.tStep;
                    
                case "vorty"
                    parameter(ii) = 3;
                    grid{ii}.FrameRate = 1/im{ii}.tStep;
                otherwise
                    parameter = NaN;
                    grid{ii}.FrameRate = 1/im{ii}.tStep;
                    
            end
            grid{ii}.parameter = parameter;
        end
        
        
    else
        parameter(ii) = grid{ii}.parameter;
        x{ii} = grid{ii}.x; y{ii} = grid{ii}.y;
        X{ii} = grid{ii}.X; Y{ii} = grid{ii}.Y;
        x2{ii} = grid{ii}.x2; y2{ii} = grid{ii}.y2;
    end
end

%% Plot the background
if isempty(x{1})
    x{1} = [0 1];
    warning('X Grid not defined')
end
if isempty(y{1})
    y{1} = [0 1];
    warning('Y Grid not defined')
end
set(ax, 'NextPlot', 'replacechildren');
switch lower(BackgroundType)
    case 'none'
        
    case {'timeseries', 'image', 'vorty', 'u', 'v'}
        BG_data = im{1}.cdata(:, :, parameter(1));
        if isfield(ArgsIn, 'Ice') && ArgsIn.Ice
            % Add in ice and pycnocline boundaries
            im_ice_threshold = ArgsIn.ice_threshold; % Threshold for brightness attributed to "Ice"
            ice_thickness = ArgsIn.ice_thickness; % Max thickness of ice from surface
            II = find_boundaries(ArgsIn.IceImageFname, 'ice', false, ice_thickness, im_ice_threshold);
            y_ice = interp1(II.x, II.y_ice, X{1}(1, :));
            y_ice(isnan(y_ice)) = max(y{1});

            ice_mask = false(im{1}.nx, im{1}.ny);
            for i = 1:im{1}.nx
                ice_mask(i, Y{1}(:, i) > y_ice(i)) = 1;
            end
            BG_data(ice_mask') = NaN;
           
        end
        imagesc(ax, x{1}, y{1}, BG_data) % Plot parameter as colour plot
        try
            caxis(im{1}.caxis')
        end
        
        if strcmpi(BackgroundType, 'vorty')
            caxis(ax, [-6 6])
            colormap(ax, newbluewhitered);
            c = colorbar;
            ylabel(c, '$\omega (s^{-1})$', 'interpreter' , 'latex')
        elseif strcmpi(BackgroundType,'u') || strcmpi(BackgroundType, 'v')
            caxis(ax, [-0.15 .15])
            colormap(ax, cmocean('balance'));
            c = colorbar;
            ylabel(c, ['$', BackgroundType, ' (ms^{-1})$'], 'interpreter' , 'latex')
        else
            if isfield(ArgsIn, 'colormap')
                colormap(ax, ArgsIn.colormap)
            else
                
                colormap(ax, singlecycle);
            end
        end
    case 'speed'
        BG_data = sqrt(im{1}.cdata(:, :, 1).^2 + im{1}.cdata(:, :, 2).^2);
        imagesc(ax, x{1}, y{1}, BG_data) % Plot parameter as colour plot
        caxis(ax, [0 .15])
        colormap(ax, cmocean('speed'));
        c = colorbar;
        ylabel(c, ['$ |u, v| (ms^{-1})$'], 'interpreter' , 'latex')
    case 'dissipation'
        U = im{1}.cdata(:, :, 1); V = im{1}.cdata(:, :, 2);
        [~, BG_data] = dissipation_gradient_2D(U, V, grid{1}.dx, grid{1}.dy);
        
        if isfield(ArgsIn, 'Ice') && ArgsIn.Ice
            % Add in ice and pycnocline boundaries
            im_ice_threshold = ArgsIn.ice_threshold; % Threshold for brightness attributed to "Ice"
            ice_thickness = ArgsIn.ice_thickness; % Max thickness of ice from surface
            II = find_boundaries(ArgsIn.IceImageFname, 'ice', false, ice_thickness, im_ice_threshold);
            y_ice = interp1(II.x, II.y_ice, X{1}(1, :));
            y_ice(isnan(y_ice)) = max(y{1});

            ice_mask = false(im{1}.nx, im{1}.ny);
            for i = 1:im{1}.nx
                ice_mask(i, Y{1}(:, i) > y_ice(i)) = 1;
            end
            BG_data(ice_mask') = NaN;
           
        end
        imagesc(ax, x{1}, y{1}, log10(BG_data)) % Plot parameter as colour plot
        if isfield(ArgsIn, 'Ice') && ArgsIn.Ice
            hold on
            plot(X{1}(1, :), y_ice, 'k-')
            hold off
        end
        colormap(ax, cmocean('amp'));
        caxis(ax, [-7 -2])
        c = colorbar;
        ylabel(c, '$log_{10}(\epsilon [m^2s^{-3}])$', 'interpreter' , 'latex')

end

hold on
set(ax,'YDir','normal');
set(ax, 'XDir', 'reverse');
if ~strcmpi(BackgroundType, 'timeseries')
    daspect(ax, [1 1 1])
end

%% Add in background
if length(fnm_in) > 1
    fnm_index = 2;
else
    fnm_index = 1;
end
switch lower(ForegroundType)
    case 'quiver'
        U = im{fnm_index}.cdata(:, :, 1);
        V = im{fnm_index}.cdata(:, :, 2);
        u = U(1:dn:size(grid{fnm_index}.X, 1), 1:dn:size(grid{fnm_index}.X, 2));
        v = V(1:dn:size(grid{fnm_index}.X, 1), 1:dn:size(grid{fnm_index}.X, 2));
        vekLeg('southwest', q_scale, .1, 'k');
        
        vekplot2(grid{fnm_index}.x2,grid{fnm_index}.y2,u,v,q_scale, 'k')
    case 'streamline'
        U = im{fnm_index}.cdata(:, :, 1);
        V = im{fnm_index}.cdata(:, :, 2);
        %hh = streamline(grid{fnm_index}.X, grid{fnm_index}.Y, U, V, grid{fnm_index}.x2, grid{fnm_index}.y2);
        hh = mstreamline(grid{fnm_index}.X, grid{fnm_index}.Y, U, V, dn);
        % TODO: Fix this mstreamline function properly!
        set(hh,'Color','k')
        
end
if isfield(ArgsIn, 'x1')
    hold on
    area(ArgsIn.x1, ArgsIn.y1, -.3 ,'FaceColor', 'w');
end

xlim([min(grid{1}.x) max(grid{1}.x)])
ylim([min(grid{1}.y) max(grid{1}.y)])
set(ax,'YDir','normal');
set(ax, 'XDir', 'reverse');
if ~strcmpi(BackgroundType, 'timeseries')
    daspect(ax, [1 1 1])
end
hold off
%% Vorticity Stats
if nargout >  3
    Stats.MaxVorticity = max(data, [], 'all'); % Maximum
    Stats.MinVorticity = min(data, [], 'all'); % Minimum
    [Stats.MaxAbsVorticity] = max(abs(data(:))); % Maximum absolute vorticity
    Stats.MinAbsVorticity = min(abs(data(:)));
    [ind, ~] = find(data==Stats.MaxAbsVorticity | data==-Stats.MaxAbsVorticity); % x location of max abs vorticity
    Stats.LocMaxVorticity = grid.xi(ind);
    
    %% Horizontal Velocity Stats
    u_data = im.cdata(:, :, 1);
    
    Stats.MaxU = max(u_data, [], 'all'); % Maximum
    Stats.MinU = min(u_data, [], 'all'); % Minimum
    [Stats.MaxAbsU] = max(abs(u_data(:))); % Maximum absolute vorticity
    Stats.MinAbsU = min(abs(u_data(:)));
    [ind, ~] = find(data==Stats.MaxAbsU | data==-Stats.MaxAbsU); % x location of max abs vorticity
    Stats.LocMaxU = grid.xi(ind);
    
    v_data = im.cdata(:, :, 2);
    
    Stats.MaxV = max(v_data, [], 'all'); % Maximum
    Stats.MinV = min(v_data, [], 'all'); % Minimum
    [Stats.MaxAbsV] = max(abs(v_data(:))); % Maximum absolute vorticity
    Stats.MinAbsV = min(abs(v_data(:)));
    [ind, ~] = find(data==Stats.MaxAbsV | data==-Stats.MaxAbsV); % x location of max abs vorticity
    Stats.LocMaxV = grid.xi(ind);
end
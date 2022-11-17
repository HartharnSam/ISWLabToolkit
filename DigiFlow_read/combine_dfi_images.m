function combine_dfi_images(input,  output)
%COMBINE_DFI_IMAGES - Combines multiple cameras onto one axes
% Input .dfi file either PIV output, or output from Transform to WCS
% uncompressed and uncompacted
%
% Syntax:  combine_dfi_images(input, output)
%
% Inputs:
%    input - Structure containing information relevant to the input data
%    output - Structure containing information relevant to the output
%    plot/video (see example)
%
% Example:
%    clc; clearvars; close all;
%    % Set input structure 
%    input.Type = 'Experimental'; % Experimental for raw image, or "PIV"
%    input.Ax = subplot(1, 1, 1); % Optional axes input
%    input.n_cameras = 2;
%    input.fnames{1} = 'CamA/output_####.dfi';
%    input.fnames{2} = 'CamB/output_####.dfi';
%    input.placeholder = '####';
%    
%    % Set output structure
%    output.xlimits = [4.5 5.5]; % WCS x limits of output 
%    output.ylimits = [0 .3]; % WCS ylimits of output
%    output.XDir = 'reverse'; % Direction of x axis (both this and YDir are optional);
%    output.t = [1:300]; % Vector (or single number) of frame numbers to output
%    output.fname = 'test.mp4'; % Filename of output, could be an image (e.g. test.png);
% 
%    % Run the script
%    combine_dfi_images(input, output);
%
% Other m-files required: subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME2
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 09-Sep-2021; Last revision: 09-Sep-2021
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
% Sort some parameters out based on what the input will show
switch (input.Type)
    case 'Experimental'
        input.parameter = 1;
        if ~isfield(input, 'Colormap')
            input.Colormap = 'singlecycle';
        end
        input.ColorRange = [0 1]*255;
    case 'PIV'
        input.parameter = 2;
        input.Colormap = 'newbluewhitered';
        input.ColorRange = [-1 1]*6;
    otherwise % means side is empty - do nothing
end

% Open figure
fig = gcf;
if ~isfield(input, 'Ax')
    input.ax = subaxis(1, 1, 1);
end
input.Ax.XLimMode = 'manual';
xlim(input.Ax, output.xlimits);
ylim(input.Ax, output.ylimits)
if isfield(output, 'YDir')
    set(input.Ax ,'YDir', output.YDir);
end
if isfield(output, 'XDir')
    set(input.Ax , 'XDir', output.XDir);
end

% Start plotting
if length(output.t)>1
    vid = VideoWriter(output.fname, 'MPEG-4');
    vid.Quality =100;
    open(vid);
    doNextLoop = true;
end

for t = output.t
    if ~doNextLoop % Break out loop if end of files
        disp('End of files, finished video early')
        break
    end
    
    for i = 1:input.n_cameras
        framenumber = t;
        frame_fname = convertCharsToStrings(strrep(input.fnames{i}, input.placeholder, sprintf('%04d', framenumber)));
        
        % Plot background
        im = dfireadvel(frame_fname);
        if t == output.t(1) %Collect grid info
            input.Grid{i} = dfi_grid_read(im);
            %input.Grid{i} = struct('dx', dx, 'dy', dy, 'nx',nx, 'ny', ny, 'x', x, 'y', y);
            figure(fig);
        end
        imagesc(input.Ax, input.Grid{i}.x, input.Grid{i}.y, im.cdata(:, :, input.parameter), input.ColorRange);
        hold(input.Ax, 'on')
        % Check if there's a next loop
        frame_fname = convertCharsToStrings(strrep(input.fnames{i}, input.placeholder, sprintf('%04d', framenumber+1)));
        if ~isfile(frame_fname)
            doNextLoop = false;
        end
    end
    eval(['colormap(input.Ax, ', input.Colormap, ');']);
    
    box on
    xlim(input.Ax, output.xlimits);
    ylim(input.Ax, output.ylimits);
    set(input.Ax, 'YDir', 'normal');
    set(input.Ax, 'XDir', 'reverse');
    hold(input.Ax, 'off');
    xlabel('x (m)'); ylabel('y (m)');
    if isfield(output, 'daspect')
    daspect(output.daspect)
    end
    if length(output.t)>1
        writeVideo(vid, getframe(input.Ax));
        completion(t-output.t(1), output.MaxFrames-output.t(1));
    end
end

if length(output.t)>1
    close(vid)
else
    [~, ~, filetype] = fileparts(output.fname);
    print(output.fname, ['-d' filetype]);
end

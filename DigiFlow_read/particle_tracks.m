function particle_tracks(filename, max_t, isPlot)
%PARTICLE_TRACKS - reads in ptv data from each timestep incrementally, and
%collates for each particle track, saves as ptv_tracks.mat
%Can also track numbers vs locations over time if isPlot is "true" (doesn't save this)
%
% Other m-files required: ptv_read, completion
% Subfunctions: none
% MAT-files required: none
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 26-Nov-2020; Last revision: 08-Feb-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

%% Run the function
if nargin < 3
    isPlot = false;
end

max_files = numel(dir([filename, '/ptv_bas*']))-1;
if nargin<2 || max_files<max_t
    max_t = max_files;
end

for i = 1:max_t
    output = ptv_read(filename, i-1, 'full');
    emptyoutput = [];
    for j = 1:length(output)
        if ~isempty(output{j})
            emptyoutput(end+1) = j;
        end
    end
    if i == 1
        
        % Put something into the ptv structure to identify variables
        ptv.Variables = fieldnames(output{emptyoutput(1)});
        ptv.Variables{end+1} = 'x_Area';
        ptv.Variables{end+1} = 'y_Area';
        ptv.data = cell(1, 1);
    end
    
    for j = 1:length(emptyoutput)
        ind = emptyoutput(j);
        try
            if (output{ind, 1}.idTrack > length(ptv.data)) || isempty(ptv.data{output{ind, 1}.idTrack})
                ptv.data{output{ind, 1}.idTrack} = nan(max_t, length(ptv.Variables));
            end
            
            % Calculate the area centroids
            output{ind}.x_Area = output{ind}.x - output{ind}.mismatchx;
            output{ind}.y_Area = output{ind}.y - output{ind}.mismatchy;
            
            % Store each particle's data
            
            ptv.data{output{ind, 1}.idTrack}(i,:) = cell2mat(struct2cell(output{ind}));
            
        catch
           % error(['Error on processing timestep ' num2str(i)])
        end
        if isPlot
            text(output{ind, 1}.x, output{ind, 1}.y, num2str(output{ind, 1}.idTrack));
            hold on
        end
        
    end
    
    if isPlot
        drawnow
        hold off
        xlim([0 1000])
        ylim([0 1000])
        pause(.1)
        clf
    end
    completion(i, max_t);
end
ptv.n_particles = length(ptv.data);
ptv.n_timesteps = i;

save(fullfile(filename, 'ptv_tracks.mat'), 'ptv');
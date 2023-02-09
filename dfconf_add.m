function dfconf_add(paramName, paramValue)
%DFCONF_ADD - Add (or replace) parameters to the df.conf file
%
% Example:
%    dfconf_add('A_w', 0.073);
%
% Other m-files required: dfconf_params
% Subfunctions: none
% MAT-files required: none
%
% See also: dfconf_params
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 30-Jan-2023; Last revision: 30-Jan-2023
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

% First see if there exists a df.conf file
if isfile('./df.conf')
    params = dfconf_params;
else
    params =  struct();
end


% Then see if the given parameter exists already (replacement) or not
% (addition)
if isa(paramValue, 'double');
    str = sprintf('%s = %6.5f', paramName, paramValue);
else %Assume it's a character
    str = sprintf('%s = %s', paramName, paramValue);
end

if isfield(params, paramName)

    fid = fopen('df.conf');
    data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
    fclose(fid);
    % Read in file and find where paramName is found
    for I = 1:length(data{1})
        tf = strncmp(data{1}{I}, paramName, length(paramName));
        if tf == 1
            data{1}{I} = str;
        end
    end
    
    % Re-write df.conf file with substitution
    fid = fopen('df.conf', 'wt');
    for I = 1:length(data{1})
        fprintf(fid, '%s \n', char(data{1}{I}));
    end

else
    fid = fopen('df.conf', 'at');
    fprintf(fid, '%s \n', str);
end

fclose(fid);
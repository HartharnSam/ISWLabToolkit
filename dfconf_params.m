function params = dfconf_params()
%  DFCONF_PARAMS   Parses the df.conf file into a structure.
%
%  Usage:
%    params = dfconf_params()
%
%  Outputs:
%    params	- a structure containing fields and values from df.conf
%
%  Sam Hartharn-Evans, from spins_params (David Deepwell, 2015)

try
    fileID = fopen('df.conf');
    tline = fgetl(fileID);
    while ischar(tline)
        if ~isempty(tline) && ~strcmp(tline(1),'#')
            lin = strtrim(strsplit(tline,'='));
            nam = lin{1};
            val = lin{2};
            if is_str_numeric(val)
                val = str2double(val);
            end
            params.(nam) = val;
        end
        tline = fgetl(fileID);
    end
    fclose(fileID);
    %%
catch
    warning('df.conf was not found or has a parameter with no value.')
    params = struct();
end
end

function out = is_str_numeric(s)
% boolean output whether string s is numeric
% from http://rosettacode.org/wiki/Determine_if_a_string_is_numeric#MATLAB
    out = ~isempty(parse_float(s));
end

% Returns the float (double) if true, empty array otherwise.
function f = parse_float(s)
    [f_in_cell, pos] = textscan(s, '%f');
    % Make sure there are no trailing chars. textscan(..) is greedy.
    if pos == length(s)
        f = f_in_cell{:};
    else
        f = [];
    end
end
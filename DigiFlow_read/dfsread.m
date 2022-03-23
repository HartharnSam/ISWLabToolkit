function dfsread(filename)
% Incomplete function to try and read in data from DigiFlow_Status files

if nargin == 0 
    filename = 'DigiFlow_Status.dfs';
end
%%

%% Setup the Import Options
dataLines = [1, Inf];

opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "";

% Specify column names and types
opts.VariableNames = "Main";
opts.VariableTypes = "char";
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
status = readtable(filename, opts);

%% Find wcs
CoordSystemInd = find(ismember(status.Main, {'# Coordinate Systems'}));
CoordSystemEnds = find(ismember(status.Main, {'# Set default coordinate system'}));
DefaultCoordSystem = status.Main{CoordSystemEnds+1};
DefaultCoordSystem = DefaultCoordSystem((length('coord_system_set_default(name:= "')):end-3);

if length(CoordSystemInd) == 1
    CoordSystemStartEnds = [CoordSystemInd CoordSystemEnds];
else
    CoordSystemStartEnds = [CoordSystemInd; CoordSystemInd(2:end) CoordSystemEnds];
end

for i = 1:length(CoordSystemInd)
    
    CoordSystemNametmp = status.Main{CoordSystemInd(i)+1};
    CoordSystemName{i} = CoordSystemNametmp((length('# Coordinate System: ')+1):end);
    
    
    %inds = find(ismember(status.Main, {['coord_system_add_point(name:="', CoordSystemName{i}]}));
end


    
    


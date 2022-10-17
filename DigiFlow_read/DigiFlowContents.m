function contents = DigiFlowContents(category)
%DigiFlowContents - Lists functions of the DigiFlow_dfi_read package
%
% Inputs:
%    category - sub-folder name of the package
%
% Outputs:
%    contents - cell array of function names
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: SPINSStartup
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 21-Feb-2022; Last revision: 21-Feb-2022
% MATLAB Version: 9.10.0.1851785 (R2021a) Update 6

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%
%% Analysis
% if nargin == 0
%     open DigiFlowContents
%     return
% end
if nargin == 0
    category = '';
end
m_path = mpath;
string_chars = ' \n ';
content = ls([m_path, category, '\*.m']);
fprintf(['\n ------------------------ \n ' category ' \n ------------------------ \n'])
for i = 1:size(content, 1)
    part = content(i, :);
    [~, part] = fileparts(part);
    string_chars = [string_chars part]; %#ok<AGROW> 
    if mod(i-1, 4) == 1
        string_chars = [string_chars ' \n ']; %#ok<AGROW> 
    else
        string_chars = [string_chars ' \t ']; %#ok<AGROW> 
    end
end
string_chars = [string_chars ' \n '];

fprintf(string_chars);

if nargout == 0
return 
else
    contents = content;
end
end

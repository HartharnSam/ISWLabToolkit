function m_path = mpath
% Returns the full path to the function that called mpath
% Made to make constant relative file paths (relative to the function)
% a bit easier
%
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 07-Feb-2022; Last revision: 07-Feb-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)

st = dbstack('-completenames', 1);
m_path = fileparts(st(1).file);
m_path = [m_path, '/'];
end
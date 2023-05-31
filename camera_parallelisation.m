%CAMERA_PARALLELISATION - Checks cameras are parallel to tank from
%calibration measurements (detailed below)
%
% Set up - rule centred at i = 512 pixel (mid image), 
%   Lower Bound is i value of left hand end of rule, 
%   Upper Bound is i value of right hand end of rule
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 09-‎Mar-‎2020; Last revision: 19-Jan-2022
% MATLAB Version: 9.10.0.1602886 (R2021a)
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------clc; clearvars; close all
centre = 512;
margin = 1;
p1 = input('Lower Point = ');
p2 = input('Upper Point = ');

p1_c = centre - p1;
c_p2 = p2 - centre;

if p1_c - c_p2 > margin % (i.e. LHS bigger than RHS by more than theshold amount)
    disp('Turn Camera Counter-clockwise')
        disp(['Difference = ', num2str(p1_c), ' vs ', num2str(c_p2)]);

elseif c_p2 - p1_c > margin % (i.e. RHS bigger than LHS by more than threshold amount)
    disp('Turn Camera Clockwise')
        disp(['Difference = ', num2str(p1_c), ' vs ', num2str(c_p2)]);

else
    disp('Camera parallel')
end


% Example run script to calculate and then plot dissipationy things

clc; clearvars; close all; 
piv_dir = 'D:\OneDrive - Newcastle University\02_PhD_Project\04_Ice_Covered_Waters\02_Raw_data\01_CameraData\050122\CamB';
piv_fname = ['\piv_0590.dfi'];
image_dir = piv_dir;
image_fname = ['\output_0590.dfi'];
diagnostic = false;

ice_props.im_ice_threshold = 100; % Threshold for brightness attributed to "Ice"
win_size = 64; 
im = dfireadvel([piv_dir piv_fname]);
dx = abs(im.xWorldPerPixel); dy = abs(im.yWorldPerPixel);
sample_freq = 2*pi*dx; f_nyq = sample_freq/2; 
Lx = abs(im.nx *im.xWorldPerPixel);
k_fit_range = [2*pi/Lx 1/f_nyq];
nu =1.7E-6;
ice_props.ice_thickness = .05; % Max thickness of ice from surface

[grd, e_dir, II] = calc_dissipation(piv_dir, piv_fname, win_size, nu);
%%

ice_lines = [true, true, true, false]; % Which plots are going to have the presence of "ice" added to them?
settings.x_lim = [4.11 4.65];
settings.y_lim = [0 .3];
settings.diss_lim = [-8 -3];
settings.image_fname = image_fname; settings.piv_file = piv_fname;
settings.image_dir = image_dir; settings.piv_dir = piv_dir;

parameters = {'direct_windowed', 'spectral'};


plot_dissipation(parameters, settings, grd, e_dir, II, ice_lines);

% TODO LOTS
% Then plot?
function convert_video(input_fnm, output_fnm)
% Converts video (typically from .avi) to .mp4 format. Doesn't keep any
% audio
% INPUTS:
%   - input_fnm : Filename + path of the input video (e.g.
%   'C:\Users\Sam\testideo.avi')
%   - input_fnm : Filename + path of the output video (e.g.
%   'C:\Users\Sam\testideo.mp4')
%
% See also: 
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 04-Oct-2021; Last revision: 04-Oct-2021
% MATLAB Version: 9.10.0.1602886 (R2021a)

%%%%%%%%%%%%%%%%%%
%%%% RUN CODE %%%%
%%%%%%%%%%%%%%%%%%
clc; close all; % To make sure any videos are shut
videoReader = VideoReader(input_fnm);
videoFWriter = VideoWriter(output_fnm, 'MPEG-4');
videoFWriter.FrameRate = videoReader.FrameRate;
videoFWriter.Quality = 100;

open(videoFWriter);
max_count = [videoReader.NumFrames];
n_frames = max(videoReader.NumFrames, max_count);

%% Load and save video frame by frame
count= 1;
while hasFrame(videoReader) && count<max_count % checks both that frame less than initial "max frames" and that there is another frame to read in the video
    videoFrame = readFrame(videoReader); % Read input video frame

    writeVideo(videoFWriter, videoFrame); % Export output video frame
    completion(count, max_count, .05);
    count = count+1;
end

try
    release(videoReader);
end
close(videoFWriter);
disp('Complete')
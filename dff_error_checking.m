close all; clc; clearvars;
fid = fopen('original.dfm.dff');
i = 1;
lines = string.empty;
while ~feof(fid)
    tline = fgetl(fid);
    tline = strtrim(tline);
    if isempty(tline) || tline(1) == '%' || strcmpi(tline, '];')
        continue
    end
    line = split(tline, ':=');
    lines(i,1) = strtrim(line{1});
    lines(i,2) = strtrim(line{2});
    i = i + 1;
end

% Put the data into a MATLAB structure
for i = 1:length(lines)
    fieldname = split(lines(i, 1), '.');
    strct.(fieldname(2)) = double(split(erase(lines(i, 2), '[')));
    % Note I've butcherred this. It does nothing for any non-array data,
    % nor for leading spaces. HOWEVER it works for this purpose just fine
end

plot(strct.acqFrame, strct.acqTime);
hold on
plot(strct.acqFrameCaptured, strct.saveTime);
plot(strct.acqFrameCaptured, strct.saveTimeCaptured);

xline(strct.missFrame(~isnan(strct.missFrame)));

xlabel('Frame Number'); ylabel('Time (s)')
legend('Acq Time', 'SaveTime', 'SaveTimeCaptured', 'MissedFrames', 'Location', 'best')

fprintf('FPS Achieved %2.3f \n', strct.fpsAchieved(~isnan(strct.fpsAchieved)));
fprintf('Missed Frames %i of  total %i \n', length(find(~isnan(strct.missFrame))), length(find(~isnan(strct.acqTimeCaptured))));




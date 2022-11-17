clc; clearvars; close all; 

ice_threshold = 160/255;
figure(1);
im = dfireadvel('ice_ts.dfi');
grid = dfi_grid_read(im);
data = im.cdata(:, :, 1);
data(data<ice_threshold) = NaN;

imagesc(grid.x, grid.y, data);
colormap('gray')

%% Detect peaks
figure(3); clf
figure(2); clf;

xpeaks = nan(1, grid.ny); widths = xpeaks;
for ii = 1:grid.ny

    %figure(2);
    %plot(grid.xi, data(ii, :));
    
    % Identify the peak initially
    [~, xpeaks_tmp, widths_tmp] = findpeaks(flip(data(ii, :)), flip(grid.xi),'SortStr','descend',...
        'NPeaks', 5);
    % TODO: Add in flat peak detection, and adjust the location of the flat
    % to the centre instead?

    % Make a prediction based on the previous timestep
    if ii > 1 && ~isnan(xpeaks(ii-1)) length(xpeaks_tmp)>=1;
        region = -25*widths(ii-1).^2 + widths(ii-1)
        xpeak_min = xpeaks(ii-1)-region/2;
        xpeak_max = xpeaks(ii-1)+region/2;
        
        if xpeaks_tmp(1)<xpeak_min || xpeaks_tmp(1)>xpeak_max
            xpeakind = find(xpeaks_tmp>xpeak_min & xpeaks_tmp<xpeak_max);
            if length(xpeakind)>=1
                xpeaks_tmp(1) = xpeaks_tmp(xpeakind(1));
                widths_tmp(1) = widths_tmp(1);
            else
                warning(['On ii = ', num2str(ii), ' Prediction failed!'])
                xpeaks_tmp(1) = NaN;
                widths_tmp(1) = NaN;
            end

        end
    end
    if ~isempty(xpeaks_tmp)
        xpeaks(ii) = xpeaks_tmp(1);
        widths(ii) = widths_tmp(1);
        hold on
        plot(xpeaks(ii), .9, 'kx')
        hold off
    else
        xpeaks(ii) = NaN; widths(ii) = NaN;
    end

    %xpeaks{ii} = grid.xi(peaks{ii});


    %pause(.05);

end

figure(1);
hold on
plot(xpeaks, grid.yi, 'rx')

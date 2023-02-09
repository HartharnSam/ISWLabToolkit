function fig = plot_dissipation(parameters, settings, e_spec, e_dir, II, ice_lines)
%PLOT_DISSIPATION - Make some plots of input dissipation, calculated by
%calc_dissipation
%
% Syntax:  fig = plot_dissipation(parameters,settings,grd, e_dir, II, ice_lines)
%
% Inputs:
%    parameters - 
%    settings - structure containing plot settings (y_lim, x_lim,
%    diss_lim, image_dir, image_file)
%    e_spec - e_spec output by calc_dissipation
%    e_dir - e_dir output by calc_dissipation
%    II - Structure identifying boundaries, as produced by find_boundaries
%    ice_lines - switch of which plots to add in ice masking to
%   
% Outputs:
%    fig - Figure handle for plotted figure
%
% Other m-files required: cmocean, dfireadvel, newbluewhitered, plot_dfi, subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: RUN_DISSIPATION,  CALC_DISSIPATION
% Author: Sam Hartharn-Evans, from Peter Sutherland's scripts
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 21-Sep-2021; Last revision: 21-Oct-2022
% MATLAB Version: 9.12.0.2009381 (R2022a) Update 4

set(0,'defaultAxesFontSize',12)
set(0,'defaultAxesLineWidth',1.0)
set(0,'defaultTextFontSize',12)
set(0,'defaultLineLineWidth',1.0)

%% Set up ready for variable plots
%% Make the Plots
%% Identify which parameters to be used & their positions
n_plots = length(parameters);
n_columns = 2;
n_rows = ceil(n_plots/n_columns);
% List of possible parameters, identifies the order
positioning.direct_windowed = find(strcmpi(parameters, 'direct_windowed'));
positioning.direct_unwindowed = find(strcmpi(parameters, 'direct_unwindowed'));
positioning.spectral = find(strcmpi(parameters, 'spectral'));
positioning.profile = find(strcmpi(parameters, 'profile'));
% Accessory ones
positioning.piv = find(strcmpi(parameters, 'piv'));
positioning.rawimage = find(strcmpi(parameters, 'rawimage'));

%% Plot Dissipation estimates
diss_lim = settings.diss_lim;
x_lim = settings.x_lim; y_lim = settings.y_lim;

fig = figure(12); clf;

%% --- Plot direct estimates before windows --- %%
if ~isempty(positioning.direct_unwindowed)
    ax_pos = positioning.direct_unwindowed;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    imagesc(e_dir.x, e_dir.y, log10(e_dir.e_image)); shading flat; cmocean('amp');
    set(ax_pos, 'YDir', 'normal', 'XDir', 'reverse');
    axis equal;
    caxis(diss_lim);
    xlim(x_lim); ylim(y_lim);
    xlabel('x [m]');
    ylabel('z [m]');
    title({'Direct \epsilon Estimate', '(non-windowed)'})
    colorbar;
    hold on
    
    if ice_lines(ax_pos)
        plot(II.qpx,II.qpy,'w','HandleVisibility','off');
        plot(II.x,II.y_ice,'w','LineWidth',4,'HandleVisibility','off');
        plot(II.x,II.y_ice,'b','LineWidth',1);
    end
    plot(II.x,II.y_interface,'w-','LineWidth',2,'HandleVisibility','off');
    plot(II.x,II.y_interface,'k--','LineWidth',1);
    
end

%% --- Plot direct estimates with windows --- %%
if ~isempty(positioning.direct_windowed)
    ax_pos = positioning.direct_windowed;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    e_dir_plot = e_dir.e_filt;
    % Use contourf with line style 'none' due to pixelation of imagesc
    contourf(e_spec.x, e_spec.y, log10(e_dir_plot), diss_lim(1):.1:diss_lim(2), 'LineStyle', 'none');
    %imagesc(grd.x, grd.y, log10(e_dir_plot));
    set(s(ax_pos), 'YDir', 'normal', 'XDir', 'reverse');
    axis equal;
    caxis(diss_lim);
    cmocean('amp');
    xlim(x_lim); ylim(y_lim);
    xlabel('x [m]'); ylabel('z [m]');
    title({'Direct \epsilon Estimate', '(windowed)'})
    c = colorbar; ylabel(c, '$log_{10}(\epsilon [m^2s^{-3}])$', 'interpreter' , 'latex');
    
    shading flat;
    hold on
    if ice_lines(ax_pos)
        
        plot(II.qpx,II.qpy,'w','HandleVisibility','off');
        plot(II.x,II.y_ice,'w','LineWidth',4,'HandleVisibility','off');
        plot(II.x,II.y_ice,'b','LineWidth',1);
    end
    plot(II.x,II.y_interface,'w-','LineWidth',2,'HandleVisibility','off');
    plot(II.x,II.y_interface,'k--','LineWidth',1);
    
    %    disp(['<direct> = ' num2str(mean(e_dir_plot(ok_comp), 'omitnan')) ]);
end
%% --- Plot spectral dissipation estimate --- %%
if ~isempty(positioning.spectral)
    ax_pos = positioning.spectral;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    spec_plot = grd.1d_Nasmyth;
    
    % pcolor(linspace(x_range(1),x_range(2),n),linspace(y_range(2),y_range(1),m),log10(e_grad_im(end:-1:1,:)))
    contourf(e_spec.x,e_spec.y,log10(spec_plot),diss_lim(1):.1:diss_lim(2), 'LineStyle','none');
    set(s(ax_pos), 'YDir', 'normal', 'XDir', 'reverse'); axis equal;
    shading flat;
    colorbar;
    caxis(diss_lim);
    cmocean('amp')
    axis equal;
    xlim(x_lim); ylim(y_lim);
    xlabel('x [m]');
    ylabel('z [m]');
    c_bar = colorbar;
    title('\epsilon from spectra');
    set(get(c_bar,'ylabel'),'string','\epsilon [m^2/s^3]','VerticalAlignment','Top','Rotation',90,'FontSize',get(0,'defaultAxesFontSize'));
    hold on;
    
    if ice_lines(ax_pos)
        plot(II.qpx,II.qpy,'w','HandleVisibility','off');
        plot(II.x,II.y_ice,'w','LineWidth',4);
        plot(II.x,II.y_ice,'b','LineWidth',1);
        
    end
    plot(II.x,II.y_interface,'w-','LineWidth',2);
    plot(II.x,II.y_interface,'k--','LineWidth',1);
    %disp(['<spectral> = ' num2str(mean(spec_grd_no_ice(ok_comp), 'omitnan')) ]);
    
end

%% --- Plot mean vertical profiles --- %%
% --- spectral and direct TKE dissipation estimates ---
if ~isempty(positioning.profile)
    ax_pos = positioning.profile;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    if ice_lines(ax_pos)
        
        e_spec.e1d(e_spec.ice_mask) = NaN;
    end
    e_dir_plot = e_dir.e_filt;
    plot(log10(mean((e_spec.e1d), 2, 'omitnan')),e_spec.y,'r','LineWidth',3);
    hold on;
    plot(log10(mean(e_dir_plot, 2, 'omitnan')),e_spec.y,'k','LineWidth',3)
    xlim(diss_lim);
    grid on;
    xlabel('\epsilon [m^2/s^3]');
    ylabel('z [m]');
    ylim(y_lim);
    %legend('Spectral','Direct','location','NW');
    set(s(ax_pos), 'XDir', 'normal');
end

if ~isempty(positioning.piv)
    ax_pos = positioning.piv;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    ArgsIn.Parameter = 3;
    ArgsIn.q_scale = .2;
    ArgsIn.Resolution = 24;
    ArgsIn.PIV_Form = 'VelArrows';
    piv_dir = settings.piv_dir;
    piv_file = settings.piv_file;
    plot_dfi(fullfile(piv_dir, piv_file), 'piv',  ArgsIn, [], s(ax_pos))
    colormap(s(ax_pos), newbluewhitered)
    caxis(s(ax_pos), [-10 10])
    axis equal;
    xlim(x_lim); ylim(y_lim);
    xlabel('x [m]');
    ylabel('z [m]');
    c_bar = colorbar;
    set(get(c_bar,'ylabel'),'string','\omega [s^{-1}]','VerticalAlignment','Top','Rotation',90,'FontSize',get(0,'defaultAxesFontSize'));
    title('Velocity/Vorticity');
    hold on;
    
    if ice_lines(ax_pos)
        plot(II.qpx,II.qpy,'w','HandleVisibility','off');
        plot(II.x,II.y_ice,'w','LineWidth',4);
        plot(II.x,II.y_ice,'b','LineWidth',1);
        
    end
    plot(II.x,II.y_interface,'w-','LineWidth',2);
    plot(II.x,II.y_interface,'k--','LineWidth',1);
end

if ~isempty(positioning.rawimage)
    ax_pos = positioning.rawimage;
    s(ax_pos) = subaxis(n_rows, n_columns, ax_pos);
    image_dir = settings.image_dir;
    image_file = settings.image_file;
    im = dfireadvel(fullfile(image_dir, image_file)); % read .dfi file
    imagesc(s(ax_pos), Grid.x_range, Grid.y_range, im.cdata(:, :, 1), [0 1]);
end
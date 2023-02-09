function [S,fx,fy] = spec2_ps_nopad(data, freq_sample_x, freq_sample_y,n_fft_x,n_fft_y,win_type)
%SPEC2_PS_NOPAD - function to compute the 2-D power spectrum of a data set
% (typically an image).
% 
% [S,fx,fy] = spec2_ps_nopad(data,freq_sample_x,freq_sample_y,n_fft_x,n_fft_y,win_type)
%
% Peter Sutherland - 14/03/2011
%
% Inputs:
%    data - Description
%    freq_sample_x - Sampling frequency in x (1/dx)
%    freq_sample_y - Sampling frequency in y (1/dy)
%    n_fft_x - Width of window in number of x pixels
%    n_fft_y - Height of window in number of y pixels
%    win_type - (String) Window type:
%               'Hanning'
%               'none'
%
% Outputs:
%    output1 - Description
%    output2 - Description

% USE EVEN WINDOW SIZES (POWERS OF 2 ARE BEST)
% Check if this is the case
if n_fft_x ~= n_fft_y
    warning('Window size should be even')
end
if mod((log2(n_fft_x)), 1) % is it a power of 2?
    warning('Window size best to be power of 2')
end
% data = XT;
% sample_freq_x = 2*pi/.25;
% sample_freq_y = 2*pi/VEL.dx_regrid;
%
% n_fft_y = 1024; n_fft_x = 1024; 
%
%if nargin < 4, n_breaks = 1; end
if nargin < 6, win_type = 'none'; end

[m,n] = size(data);
% N_s_x = floor(n/(n_breaks_x+1)*2); % length of x windows
% N_s_y = floor(m/(n_breaks_y+1)*2); % length of y windows

n_win_x = floor(2*n/n_fft_x)-1; % number of windows in each direction.
n_win_y = floor(2*m/n_fft_y)-1; 


% shift the last window if within 10% of being able to fit another.
if (n_fft_x/2 - (n - (n_fft_x/2*(n_win_x+1)))) < (0.1*n_fft_x)
    shift_last_x = -(n_fft_x/2 - (n - (n_fft_x/2*(n_win_x+1))));
    n_win_x = n_win_x+1;
else
    shift_last_x = 0;
end
if (n_fft_y/2 - (m - (n_fft_y/2*(n_win_y+1)))) < (0.1*n_fft_y)
    shift_last_y = -(n_fft_y/2 - (m - (n_fft_y/2*(n_win_y+1))));
    n_win_y = n_win_y+1;
else
    shift_last_y = 0;
end

% find start and stop of each window in each direction
win_ss_x = NaN(n_win_x,2);
for i_win_x = 1:n_win_x
    ind_start_x = (i_win_x-1)*n_fft_x/2;
    win_ss_x(i_win_x,:) = [ (ind_start_x+1) (ind_start_x+n_fft_x) ];
end
if shift_last_x~=0, win_ss_x(end,:) = win_ss_x(i_win_x,:)+shift_last_x; end % shift the last window if within 10% of being able to fit another.

win_ss_y = NaN(n_win_y,2);
for i_win_y = 1:n_win_y
    ind_start_y = (i_win_y-1)*n_fft_y/2;
    win_ss_y(i_win_y,:) = [ (ind_start_y+1) (ind_start_y+n_fft_y) ];
end
if shift_last_y~=0, win_ss_y(end,:) = win_ss_y(i_win_y,:)+shift_last_y; end % shift the last window if within 10% of being able to fit another.

% set up window function
if strcmp(win_type,'hanning') % weights based on distance from centre of window
    win_2D = 0.5*(1-cos(2*pi*(0:n_fft_y-1)'/(n_fft_y-1)))*0.5*(1-cos(2*pi*(0:n_fft_x-1)/(n_fft_x-1)));
    win_correct_2D = (8/3).^2; % 1./mean(win_hanning(:).^2); % for Hanning?
elseif strcmp(win_type,'none')
    win_2D = ones(n_fft_y,n_fft_x);
    win_correct_2D = 1; % 1./mean(win_2D(:).^2); % for Hanning?
end

%total_windows = n_breaks_x*n_breaks_y;

%data_w = NaN(N_s_y,N_s_x,total_windows);

dx = 1/freq_sample_x;
dy = 1/freq_sample_y;

x_win = n_fft_x*dx;
y_win = n_fft_y*dy;

S = zeros(n_fft_y,n_fft_x);

for ix = 1:n_win_x
    for iy = 1:n_win_y
        DATA_raw = data(win_ss_y(iy,1):win_ss_y(iy,2),win_ss_x(ix,1):win_ss_x(ix,2)); % Re-do data into windows
        
        DATA = DATA_raw - mean(DATA_raw(:)); % Remove mean flow
        
        DATA_win = DATA.*win_2D; % 
        
        FFT_ETA = fftshift(fft2(DATA_win));
        
        S_us_temp = conj(FFT_ETA).*FFT_ETA;
        
        S = S + win_correct_2D/(x_win*y_win)*(dx*dy)^2*S_us_temp;
        
    end
end

S = S/(n_win_x*n_win_y);

fx = freq_sample_x/n_fft_x*( -floor(n_fft_x/2):ceil(n_fft_x/2-1) );
fy = freq_sample_y/n_fft_y*( -floor(n_fft_y/2):ceil(n_fft_y/2-1) );

% % % % %         elseif strcmp(mre_type,'plane')
% % % % %             
% % % % %             
% % % % %             Plane.X0 = mean([X_win(:),Y_win(:),data_temp(:)]',2);
% % % % %             A_A = [(X_win(:)-Plane.X0(1)) (Y_win(:)-Plane.X0(2)) (data_temp(:)-Plane.X0(3))];
% % % % %             [A_svd_U,A_svd_S,A_svd_V] = svd(A_A,0);
% % % % %             [A_svd_s,A_svd_i] = min(diag(A_svd_S));
% % % % %             Plane.a = A_svd_V(:,A_svd_i); % Direction cosines of normal
% % % % %             % d = U(:,svd_i)*svd_s;
% % % % %             Plane.P_fit = [Plane.a ; -(Plane.a(1)*Plane.X0(1) + Plane.a(2)*Plane.X0(2) + Plane.a(3)*Plane.X0(3))];
% % % % %             
% % % % %             mean_remove = -(Plane.P_fit(1)*X_win + Plane.P_fit(2)*Y_win + Plane.P_fit(4))/Plane.P_fit(3);
% % % % %             
% % % % %             %             figure(3); clf;
% % % % %             %             plot3(X_win(:),Y_win(:),data_temp(:),'.');
% % % % %             %             hold on;
% % % % %             %             plot3(Plane.X0(1),Plane.X0(2),Plane.X0(3),'or','MarkerSize',15);
% % % % %             %             plot3(X_win(:),Y_win(:),mean_remove,'r+');
% % % % %             
% % % % %         elseif strcmp(mre_type,'quadratic')
% % % % %             % Removes a quadratic fit to the data...
% % % % %             V = [X_win(:).^2 Y_win(:).^2 X_win(:) Y_win(:) ones(length(X_win(:)),1)];
% % % % %             abc = V\data_temp(:);
% % % % %             mean_remove = abc(1)*X_win.^2 + abc(2).*Y_win.^2 + abc(3)*X_win + abc(4)*Y_win + abc(5);
% % % % %             
% % % % %         end



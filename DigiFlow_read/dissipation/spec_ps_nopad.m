function [freq,sx,err_low_hi,dof] = spec_ps_nopad(data,sample_freq,n_breaks,win_type,conf_level)
% spec_ps.m - produces a power spectra by cutting the data up into small
% sections, windowing them, and taking the ffts of those sections. 
%
% USAGE:
%   [freq,sx,err_low_hi] = spec_ps_nopad(data,sample_freq,n_breaks,win_type,conf_level)
%
% EXAMPLE:
%   [f,St] = spec_ps_nopad(data,sample_freq,11,'hanning');
%
% sample_freq = frequency at which data was sampled.
%
% data = the time series to be studied; can be a matrix, but
%        DATA MUST BE IN COLUMNS
%
% n_fft = the number of points each fft computes
%         (powers of 2 will make it run faster)
%
% n_breaks = number of overlapping sections the data is divided
% into before windowing (must be odd)
%
% sx*err_low_hi(1) <= Spec true <= sx*err_low_hi(2)
% with a confidence level of << conf_level*100 >> percent
% (see Power spectral estimation, National semiconductor application note
% 255, November 1980)
%

% Peter Sutherland
% 24 Jun 2015 - corrected error estimation
% 23 Sep 2013 - now works with data in matrix form
% 16 Nov 2009 - now works with data in matrix form
% 10 Jun 2009 - forked from psd_si
%  3 Mar 2008 - Produces correct energy density, but normalisation scheme
%               is still a bit of a hack
% 27 Jan 2008 - Generally Works - NEEDS PROPER NORMALISATION
% 18 Dec 2007 - NOT WORKING YET

[m n] = size(data);
f_nyq = sample_freq/2;

if m==1&&n>1
    data = data';
    [m n] = size(data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3, n_breaks = 1; end

if nargin < 4, win_type = 'none'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



data = detrend(data); %-ones(m,1)*mean(data,1); % removes mean

n_s = floor(m/(n_breaks+1)*2); % length of windows


% bin contains start and stop of each window
bin = NaN(n_breaks,2);
for k = 1:(n_breaks+1)/2
    bin(k,:) = [(k-1)*n_s+1 k*n_s];
end
for k = 1:(n_breaks-1)/2
    bin((n_breaks+1)/2+k,:) = [(k-1/2)*n_s+1 (k+1/2)*n_s];
end
bin(bin<1) = 1; bin(bin>m)=m; bin = round(bin); % rounds bin indices to nearest whole number
bin = sortrows(bin,1);                          % sort bin rows

data_w = NaN(n_s,n,n_breaks);
var_dat = NaN(1,n,n_breaks);
for k = 1:n_breaks
    dw_temp = detrend(data(bin(k,1):bin(k,2),:));
    var_dat(1,:,k) = var(dw_temp,0,1);
    if strcmp(win_type,'none')
        % *
        data_w(:,:,k) = dw_temp;
    else
        % *
        data_w(:,:,k) = dw_temp.*((window(win_type,n_s)*ones(1,n)));
    end
end

n_fft = size(data_w,1);
freq = (linspace(0,f_nyq,ceil(n_fft/2)))'; % frequency vector

r_fft = fft(data_w,[],1);

sx_uf = conj(r_fft).*r_fft; % /(m^2);

sx_uf_nor = NaN(size(sx_uf));

% *
f_full = linspace(0,2*f_nyq,n_fft)';
for iwin = 1:n_breaks
    % sx_uf_nor(:,:,iwin) = sx_uf(:,:,iwin).*repmat(var(data_w(:,:,iwin),0,1)./trapz(f_full,real(sx_uf(:,:,iwin)),1),n_fft,1);
    sx_uf_nor(:,:,iwin) = sx_uf(:,:,iwin).*repmat(var_dat(1,:,iwin)./trapz(f_full,real(sx_uf(:,:,iwin)),1),n_fft,1);
end
sx = 2*mean(sx_uf_nor(1:ceil(n_fft/2),:,:),3, 'omitnan');
% *

% sx_uf1 = sum(sx_uf(1:ceil(n_fft/2),:,:),3)/n_breaks;
% 
% % ensure that integral of energy = variance of signal
% % sx = sx_uf1/mean(diff(freq)).*(ones(ceil(n_fft/2),1)*(var(data,0,1)./sum(sx_uf1,1)));
% sx = sx_uf1.*(ones(ceil(n_fft/2),1)*(var(data,0,1)./trapz(freq,real(sx_uf1),1)));


if nargout>=3
    if nargin<5, conf_level = 0.95; end
    % n_breaks = degrees of freedom
    err_low_hi = [2*n_breaks/chi2inv((1-(1-conf_level)/2),2*n_breaks), 2*n_breaks/chi2inv((1-conf_level)/2,2*n_breaks)];
end

if nargout>=4
    
    N_ts = m; % number of data points in the time series.
    M_win_2 = n_s;% /2; % half-width of the window in the time domain.
    
    % c.f. Thomson and Emery, Data analysis methods in physical
    % oceanography, p479
    dof = 8/3*(N_ts/M_win_2);
    
    % C.f. earle 1996, NOAA, p7-8
    dof = 2*n_breaks/(1+ 0.4*(n_breaks-1)/n_breaks);
    
    
    % Truncated periodogram = (N/M)
    % Daniell 2*(N/M)
    % Parzen 3.708615*(N/M)
    % bartlett 3*(N/M)
    % 'hann' = 8/3*(N/M)
    % 'hamming' = 2.5164*(N/M)
    
end



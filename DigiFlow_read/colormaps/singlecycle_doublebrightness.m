function map = singlecycle_doubleBrightness(m)
%SINGLECYCLE - DigiFlow "single cycle - double brightness" color map
%
% Syntax:  singlecycle_doubleBrightness returns a colormap with the same number of colours
%          as the current figures's colormap. If no figure exists same number of
%          colors as default colormap.
%          singlecycle_doubleBrightness(m) - returns an M-by-3 matrix containing a colormap
%
% Inputs:
%    m - Length of desired colormap
%
% Outputs:
%    map - "singlecycle_doubleBrightness" colormap
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: HSV, PARULA
% Author: Sam Hartharn-Evans
% email address: s.hartharn-evans2@ncl.ac.uk
% GitHub: github.com/HartharnSam/DigiFlow_dfi_read
% April 2020; Last revision: 09/04/2020
%
% MATLAB Version: 9.8.0.1323502 (2020a)

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
if nargin < 1
    f = get(groot, 'CurrentFigure');
    if isempty(f)
        m = size(get(groot, 'DefaultFigureColormap'), 1);
    else
        m = size(f.Colormap, 1);
    end
end

values = [
    0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000;
    0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9843 0.9686 0.9529 0.9373 0.9216 0.9059 0.8902 0.8745 0.8588 0.8431 0.8275 0.8118 0.7961 0.7804 0.7647 0.7490 0.7333 0.7176 0.7020 0.6863 0.6706 0.6549 0.6392 0.6235 0.6078 0.5922 0.5765 0.5608 0.5451 0.5294 0.5137 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5098 0.5176 0.5255 0.5333 0.5412 0.5490 0.5569 0.5647 0.5725 0.5804 0.5882 0.5961 0.6039 0.6118 0.6196 0.6275 0.6353 0.6431 0.6510 0.6588 0.6667 0.6745 0.6824 0.6902 0.6980 0.7059 0.7137 0.7216 0.7294 0.7373 0.7451 0.7529 0.7608 0.7686 0.7765 0.7843 0.7922 0.8000 0.8078 0.8157 0.8235 0.8314 0.8392 0.8471 0.8549 0.8627 0.8706 0.8784 0.8863 0.8941 0.9020 0.9098 0.9176 0.9255 0.9333 0.9412 0.9490 0.9569 0.9647 0.9725 0.9804 0.9882 0.9961;
    0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5020 0.5176 0.5333 0.5490 0.5647 0.5804 0.5961 0.6118 0.6275 0.6431 0.6588 0.6745 0.6902 0.7059 0.7216 0.7373 0.7529 0.7686 0.7843 0.8000 0.8157 0.8314 0.8471 0.8627 0.8784 0.8941 0.9098 0.9255 0.9412 0.9569 0.9725 0.9882 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9922 0.9843 0.9765 0.9686 0.9608 0.9529 0.9451 0.9373 0.9294 0.9216 0.9137 0.9059 0.8980 0.8902 0.8824 0.8745 0.8667 0.8588 0.8510 0.8431 0.8353 0.8275 0.8196 0.8118 0.8039 0.7961 0.7882 0.7804 0.7725 0.7647 0.7569 0.7529 0.7608 0.7686 0.7765 0.7843 0.7922 0.8000 0.8078 0.8157 0.8235 0.8314 0.8392 0.8471 0.8549 0.8627 0.8706 0.8784 0.8863 0.8941 0.9020 0.9098 0.9176 0.9255 0.9333 0.9412 0.9490 0.9569 0.9647 0.9725 0.9804 0.9882 0.9961
    ]';

P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,m), 'linear');
%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------
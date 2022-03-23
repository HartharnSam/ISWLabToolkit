function [x, y, dx, dy, nx, ny, x0, y0, xi, yi, X, Y]  = dfi_grid_read(im)
%DFI_GRID_READ - Reads grid from im (output of dfireadvel) and provides
%data on x and y
%
% INPUTS
% im - image structure output from dfireadvel
%
% OUTPUTS
% [x y] - Corner points of x and y grids
% [dx dy] - WCS change per pixel
% [nx ny] - Number of pixels in x and y
% [x0 y0] - Origin for x and y
% [xi yi] - Full vector of x and y coordinates
% [X Y]   - X and Y grids, as by meshgrid
% If only one output called, the first output is [Grid] - Structure
% containing all of the above outputs
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Sam Hartharn-Evans (from past versions of dfi_read package)
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 08-Sept-2021; Last revision: 09-Feb-2022
% MATLAB Version: '9.4.0.813654 (R2018a)'

% Read in from im structure
x0=im.xOriginWorld;y0=im.yOriginWorld;
dx=im.xWorldPerPixel;dy=im.yWorldPerPixel;
nx=im.nx;ny=im.ny;

% Calculate two pairs of grids
x=[x0 x0+dx*nx];
y=[y0+dy*ny y0];
if nargout > 8 || nargout == 1
    xi = linspace(x(1), x(2), nx); 
    yi = linspace(y(1), y(2), ny);
end
if nargout > 10 || nargout == 1 % Make meshgrid type x and y
    try 
        X = im.x;
        Y = im.y;
    catch
        X = linspace(x0, (x0 + dx*(nx-1)), nx);
        Y = linspace(y0, (y0 + dy*(ny-1)), ny)';
        X = repmat(X, ny,1);
        Y = repmat(flipud(Y(:)),1,nx);

    end
end
if nargout == 1
    x = struct('dx', dx, 'dy', dy, 'nx',nx, 'ny', ny, 'x', x, 'y', y, 'X', X, 'Y', Y, 'xi', xi, 'yi', yi);
end
% End of function
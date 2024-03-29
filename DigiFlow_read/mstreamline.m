function [h]=mstreamline(x,y,u,v,ds)
% MSTREAMLINE - plot streamlines from velocity field
%
% h=mstreamline(x,y,u,v); attempts to plot the streamlines of a flow
% by using the boundaries as starting points.
%
% The function calls STREAMLINE using the boundaries of the 
% velocityfield as STARTX and STARTY values.
%
% h=mstreamline(x,y,u,v,ds); downsamples the START-vectors by a
% factor ds to display fewer streamlines.
%
% use set(h,'Color','blue') to change the color of the
% streamlines in the figure.
% 
% See also: STREAMLINE
%
% Copyright J. Kristian Sveen, jks@math.uio.no
% Time Stamp: 9. March 2002, 17:51
% For use with MatPIV 1.6 and subsequent versions
 
  
% First we locate the boundaries

if nargin==4
  ds=1;
end

in=find(~isnan(u));
yy = NaN(2, size(u, 1));
xx = yy;
for ii=1:size(u,1)
  ix = find(~isnan(u(ii, :)) & u(ii, :) ~= 0);
  if ~isempty(ix)
    yy(1:2, ii)=[y(ii, ix(1)+10) y(ii, ix(end)-10)];
    xx(1:2, ii)=[x(ii, ix(1)+10) x(ii, ix(end)-10)];
  else
    yy(1:2, ii)=[nan nan];
    xx(1:2, ii)=[nan nan];
  end
end
yy=[yy(:); y(1,:).'; y(:,1); y(:,end); y(end,:).'];
xx=[xx(:); x(1,:).'; x(:,1); x(:,end); x(end,:).'];
xx=xx(:); yy=yy(:);
xx=xx(1:ds:end); yy=yy(1:ds:end);

% plot the streamlines
h=streamline(x,y,u,v,xx,yy);
set(h,'Color','red');


% clear output if not called.
if nargout==0
  clear
end

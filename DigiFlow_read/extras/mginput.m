function [out1,out2,out3] = mginput(arg1)
%MGINPUT Graphical input from mouse.
%   [X,Y] = MGINPUT(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse (or by using the Arrow Keys on some
%   systems).  Data points are entered by pressing a mouse button
%   or any key on the keyboard except carriage return, which terminates
%   the input before N points are entered.
%
%   [X,Y] = MGINPUT gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = MGINPUT(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   See also ginput
%
%   Copyright (c) 1984-97 by The MathWorks, Inc.
%   $Revision: 5.19 $  $Date: 1997/04/08 06:55:47 $
%   adapted mginput version from MatPIV 1.7
%   modified pointer shape Sam Hartharn-Evans : 17-Feb-2022 : MATLAB '9.10.0.1602886 (R2021a)'

out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC') & ~strcmp(c(1:2),'MA')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') & ~strcmp(tp,'x') & ~strcmp(tp,'micro'),
    if nargout == 1,
        if nargin == 1,
            eval('out1 = trmginput(arg1);');
        else
            eval('out1 = trmginput;');
        end
    elseif nargout == 2 | nargout == 0,
        if nargin == 1,
            eval('[out1,out2] = trmginput(arg1);');
        else
            eval('[out1,out2] = trmginput;');
        end
        if  nargout == 0
            out1 = [ out1 out2 ];
        end
    elseif nargout == 3,
        if nargin == 1,
            eval('[out1,out2,out3] = trmginput(arg1);');
        else
            eval('[out1,out2,out3] = trmginput;');
        end
    end
else
    
    fig = gcf;
    figure(gcf);
    
    if nargin == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  isstr(how_many) ...
                | size(how_many,1) ~= 1 | size(how_many,2) ~= 1 ...
                | ~(fix(how_many) == how_many) ...
                | how_many < 0
            error('Requires a positive integer.')
        end
        if how_many == 0
            ptr_fig = 0;
            while(ptr_fig ~= fig)
                ptr_fig = get(0,'PointerWindow');
            end
            scrn_pt = get(0,'PointerLocation');
            loc = get(fig,'Position');
            pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
            out1 = pt(1); y = pt(2);
        elseif how_many < 0
            error('Argument must be a positive integer.')
        end
    end
    pointer = get(gcf,'pointer');
    set(gcf,'pointer','crosshair');
    fig_units = get(fig,'units');
    char = 0;
    
    while how_many ~= 0
        % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
        eval('keydown = wfbp;', 'waserr = 1;');
        if(waserr == 1)
            if(ishandle(fig))
                set(fig,'pointer',pointer,'units',fig_units);
                error('Interrupted');
            else
                error('Interrupted by figure deletion');
            end
        end
        tel=1;
        ptr_fig = get(0,'CurrentFigure');
        if(ptr_fig == fig)
            if keydown
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
                scrn_pt = get(0, 'PointerLocation');
                set(fig,'units','pixels')
                loc = get(fig, 'Position');
                pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
                set(fig,'CurrentPoint',pt);
            else
                button = get(fig, 'SelectionType');
                if strcmp(button,'open')
                    button = b(length(b));
                elseif strcmp(button,'normal')
                    button = 1;
                    pt = get(gca, 'CurrentPoint');
                elseif strcmp(button,'extend')
                    button = 2;
                    pt = get(gca, 'CurrentPoint');
                    xa=get(gca,'XLim');xa=xa(2)-xa(1);
                    ya=get(gca,'YLim');ya=ya(2)-ya(1);
                    axis([max([0,pt(1,1)-25]) min([xa,pt(1,1)+25]),...
                        max([0,pt(1,2)-25]) min([ya,pt(1,2)+25])])
                    if tel==1
                        xa2=xa; ya2=ya;
                    end
                    tel=tel+1;
                elseif strcmp(button,'alt')
                    button = 3;
                    pt = get(gca, 'CurrentPoint');
                    if tel~=1
                        xa=get(gca,'XLim');xa=xa(2)-xa(1);
                        ya=get(gca,'YLim');ya=ya(2)-ya(1);
                        axis([0.5 xa+0.5 0.5 ya+0.5]);
                        % axis([max([0,pt(1,q)-25]) min([xa,pt(1,q)+25]),...
                        %max([0,pt(1,2)-25]) min([ya,pt(1,2)+25])])
                    end
                else
                    error('Invalid mouse selection.')
                end
            end
            pt = get(gca, 'CurrentPoint');
            
            how_many = how_many - 1;
            
            if(char == 13) % & how_many ~= 0)
                % if the return key was pressed, char will == 13,
                % and that's our signal to break out of here whether
                % or not we have collected all the requested data
                % points.
                % If this was an early breakout, don't include
                % the <Return> key info in the return arrays.
                % We will no longer count it if it's the last input.
                break;
            end
            
            out1 = [out1;pt(1,1)];
            y = [y;pt(1,2)];
            b = [b;button];
        end
    end
    
    set(fig,'pointer',pointer,'units',fig_units);
    
    if nargout > 1
        out2 = y;
        if nargout > 2
            out3 = b;
            % $$$ 	if out3==2
            % $$$ 	  xa=get(gca,'XLim');xa=xa(2)-xa(1);
            % $$$ 	  ya=get(gca,'YLim');ya=ya(2)-ya(1);
            % $$$ 	  axis([max([0,out1-25]) min([xa,out1+25]),...
            % $$$ 		max([0,out2-25]) min([ya,out2+25])
            % $$$ 	elseif out3==3
            % $$$ 	  xa=get(gca,'XLim');xa=xa(2)-xa(1);
            % $$$ 	  ya=get(gca,'YLim');ya=ya(2)-ya(1);
            % $$$ 	  axis([max([0,out1-25]) min([xa,out1+25]),...
            % $$$ 		max([0,out2-25]) min([ya,out2+25])
            % $$$ 	else
            % $$$ 		axis tight
            % $$$ 		end
        end
    else
        out1 = [out1 y];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

% Remove figure button functions
fprops = {'windowbuttonupfcn','buttondownfcn', ...
    'windowbuttondownfcn','windowbuttonmotionfcn'};
fig = gcf;
fvals = get(fig,fprops);
set(fig,fprops,{'','','',''})

% Remove all other buttondown functions
ax = findobj(fig,'type','axes');
if isempty(ax)
    ch = {};
else
    ch = get(ax,{'Children'});
end
for i=1:length(ch),
    ch{i} = ch{i}(:)';
end
h = [ax(:)',ch{:}];
vals = get(h,{'buttondownfcn'});
mt = repmat({''},size(vals));
set(h,{'buttondownfcn'},mt);

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
eval(['if nargout==0,', ...
    '   waitforbuttonpress,', ...
    'else,', ...
    '   keydown = waitforbuttonpress;',...
    'end' ], 'waserr = 1;');

% Put everything back
if(ishandle(fig))
    set(fig,fprops,fvals)
    set(h,{'buttondownfcn'},vals)
end

if(waserr == 1)
    error('Interrupted');
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

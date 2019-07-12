function hh = plot2org(x,y,xy0,linespec) ;
% PLOT2ORG - scatter plot with connecting lines to origin
%    PLOT2ORG(X,Y) creates a scatter plot of Y vs. X in which each point is
%    connected to the origin (0,0). X and Y should be equal length vectors.
%    PLOT2ORG(X,Y,[X0 Y0]) in which X0 and Y0 are scalars connects each
%    point to the origin specified by (X0,Y0).
%
%    PLOT2ORG(...,'LineSpec') uses the color, marker and linestyle
%    properties specified by the string LineSpec. See PLOT for possibilities.
%
%    H = PLOT2ORG(...) returns a two-element vector with a handle to the
%    scatter points and a handle to the connecting lines.
%    
%    Example:
%      % the data
%      x = 10*rand(10,1)-5 ;
%      y = 10*rand(10,1)-5 ;
%      plot2org(x,y,'go-') ;
%      hold on ;
%      plot2org(x,y,[mean(x) mean(y)],'r.:')
%
%    See also PLOT, SCATTER, ERRORBAR
%    and XYREFLINE (Matlab File Exchange) 

% for Matlab R13
% version 1.1 (2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created: jun 2006
% 1.1 (jun 2006) Solved problem of unwanted lines when no marker symbol was
% improved (thanks to John D for pointing this out). 

error(nargchk(2,4,nargin)) ;

% set default arguments
if nargin==2,
    xy0 = [0 0] ;
    linespec = 'o-' ;
elseif nargin==3,
    if isstr(xy0),
        linespec = xy0 ;
        xy0 = [0 0] ;
    else
        linespec = 'o-' ;
    end
end

% check arguments
N = numel(x) ;
if numel(y) ~= N,
    error('plot2org:NumberOfElementsMismatch','X and Y vectors should have the same length') ;
end
if numel(xy0) ~= 2 || ~isnumeric(xy0),
    error('plot2org:WrongOriginFormat','Origin should be a 2 element vector') ;
end

% get the line specifier
[ls,col,mark,msg] = colstyle(linespec) ;
if ~isempty(msg), error('plot2org:SymbolError',msg); end
linespecP = [mark col] ;    % Use marker only on scatter part
linespecL = [ls col] ;      % Use lines only on connections to origin

% data should be a vector
x = x(:) ;
y = y(:) ;
% calculate line positions (see e.g. errorbar for similar usage)
xx0 = [x repmat([xy0(1) nan],N,1)].' ;
yy0 = [y repmat([xy0(2) nan],N,1)].' ;

% plot the points and lines
h = plot(...
    xx0(:),yy0(:),linespecL, ...
    x,y,linespecP) ;

% make sure the markers have no lines
set(h(2),'linestyle','none')
% and the connecting lines have no markers
set(h(1),'marker','none')
% and both have the same color
set(h(2),'color',get(h(1),'color')) ;

% return the handles if requested
if nargout > 0, hh = h([2 1]) ; end
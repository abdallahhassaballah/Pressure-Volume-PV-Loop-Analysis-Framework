function h = PlotPolarLambda(L,si,varargin)
% Plot theta against lambda coordinates in the 2D polar axis.
%
%   h = L.PlotPolarLambda(si,'opt1',val1,'opt2',val2,...);
%
% Input:  - si is the surface indices, or 'epi' or 'endo'.
% Options:
%   - 'color', color plot
%   - 'noplot', true | false.
%     If this is false, then there is no plotting done. The output contains
%     [X,Y] coordinates. Default is true.
%   - Other MATLAB's plot options.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% get surface
if( ischar(si) )
    if( strcmpi(si,'endo') ), si = 1;
    elseif( strcmpi(si,'epi') ), si = L.nSurfaces;
    else error('Unknown surface %s',si); 
    end
else
    if( si<1 || si>L.nSurfaces ), error('Surface index is out of range.'); end
end

% options
opt.color = 'auto';
opt.noplot = false;
plotopts = {};

% get optional arguments
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else plotopts = [plotopts varargin(i:i+1)]; end
end

% auto color
if( strcmp(opt.color,'auto') )
    if( si==1 ), opt.color = 'r';
    elseif( si==L.nSurfaces ), opt.color = 'b';
    else opt.color = 'g';
    end
end

% get indices
npts = size(L.surfacePSPoints,1)/L.nSurfaces;
idx = L.map + (si-1)*npts;

% get theta and flip around pi axes
thetas = pi - L.surfacePSPoints(idx(:),1);
lambdas = L.surfacePSPoints(idx(:),3);
[x,y] = pol2cart(thetas,lambdas);

xs = reshape(x,size(idx));
ys = reshape(y,size(idx));
if( opt.noplot )
    h.XData = xs;
    h.YData = ys;
else
    h1 = plot(xs',ys',plotopts{:},'LineStyle','-','Color',opt.color);
    hold on;
    h2 = plot(xs,ys,plotopts{:},'LineStyle','-','Color',opt.color);
    hold off;
    
    h = [h1; h2];
end
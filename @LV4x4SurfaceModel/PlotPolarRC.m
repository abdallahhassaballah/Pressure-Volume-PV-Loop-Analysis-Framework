function h = PlotPolarRC(L,si,varargin)
% Plot Y-Z in polar coordinate, so you're looking from above.
%
%   h = L.PlotPolarRC(si,'opt1',val1,'opt2',val2,...);
%
% Input:  - si is the surface indices, or 'epi' or 'endo'.
% Options:
%   - 'color', color plot
%   - 'noplot', true | false.
%     If this is false, then there is no plotting done. The output contains
%     [X,Y] coordinates. Default is true.
%   - 'scale', number
%     Apply scale, which affects the radius. Default is 1.0
%   - 'revert_y', true | false. Default is false.
%     Reverting y axis aligns with the heart angle.
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
opt.scale = 1.0;
opt.revert_y = false;
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
npts = size(L.surfacePoints,1)/L.nSurfaces;
idx = L.map + (si-1)*npts;

% get theta and flip around pi axes
%thetas = pi - L.surfacePSPoints(idx(:),1);
%lambdas = sqrt(L.surfacePoints(idx(:),2).^2 + L.surfacePoints(idx(:),3).^2);
%[x,y] = pol2cart(thetas,lambdas);

if( opt.revert_y ), rv = -1; else rv = 1; end

xs = rv * opt.scale * reshape(L.surfacePoints(idx(:),2),size(idx));
ys = opt.scale * reshape(L.surfacePoints(idx(:),3),size(idx));

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
function h = PlotSurfacePoints(L,varargin)
% Plot surface points
%
%   h = L.PlotSurfacePoints
%   h = L.PlotSurfacePoints('opt1',val1,'opt2',val2,...)
%
% Default values for sampling point markers are dots (size 5) with colors
% from red (epi) to blue (endo).
%
% Available options:
%   - 'markers', marker_type(s)
%     Marker type can be 1 or L.nSurfaces types (see doc LineSpec).
%   - 'colors', color(s)
%     Color can be 3x1 rgb color or 3xL.nSurfaces rgb colors.
%   - 'size', marker_size
%     Only one marker size that can be provided.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% default options
opt.markers = '.';
opt.colors = [...
    linspace(1,0,L.nSurfaces); ...
    zeros(1,L.nSurfaces); ...
    linspace(0,1,L.nSurfaces)];
opt.size = 5;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ) opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option %s',varargin{i}); end    
end

% adjust options
if( numel(opt.markers)==1 ), opt.markers = repmat({opt.markers},1,L.nSurfaces); end
if( size(opt.colors,2)==1 )
    opt.colors = repmat(opt.colors,1,L.nSurfaces);
end

% then plot
hs = ishold;
for i=1:L.nSurfaces
    P = L.GetSurfacePoints(i);
    h(i) = plot3(P(:,1),P(:,2),P(:,3),'Color',opt.colors(:,i),'Marker',opt.markers{i},...
        'MarkerSize',opt.size,'LineStyle','none');
    hold on;
end
if( ~hs ) hold off; end
function h = plot(L,varargin)
% Plot the LV model in 3D axis.
%
%   h = plot(L,'opt1',val1,'opt2',val2,...)
%
% This plot combines all plots from PlotSurfacePoints, PlotSurface,
% PlotWireframe and/or PlotBoundingBox.
%
% Available options:
%   - 'surface_opts', cell_array (see help PlotSurface).
%   - 'with_points', true | false. Default is true.
%   - 'point_opts', cell_array (see help PlotSurfacePoints).
%   - 'with_wireframe', true | false. Default is true.
%   - 'wireframe_opts', cell_array (see help PlotWireframe).
%   - 'with_boundingbox', true | false. Default is false.
%   - 'boundingbox_opts', cell_array (see help PlotBoundingBox).
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% default options
opt.surface_opts = {};
opt.with_points = false;
opt.point_opts = {};
opt.wireframe_opts = {};
opt.with_wireframe = true;
opt.with_boundingbox = false;
opt.boundingbox_opts = {};

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option %s.',varargin{i}); end
end

% plot surface
h.surf = L.PlotSurface(opt.surface_opts{:});
hold on;

% % plot wireframe
 if( opt.with_wireframe )
     h.wire = L.PlotWireframe(opt.wireframe_opts{:});
 end
% 
% % plot points
% if( opt.with_points )
%     h.points = L.PlotSurfacePoints(opt.point_opts{:});
% end
% 
% % plot bounding box
% if( opt.with_boundingbox )
%     h.bb = L.PlotBoundingBox(opt.boundingbox_opts{:});
% end

axis equal;
xlim(xlim + [-5 5]);
ylim(ylim + [-5 5]);
zlim(zlim + [-5 5]);
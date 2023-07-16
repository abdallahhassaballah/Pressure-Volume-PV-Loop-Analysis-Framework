function h = PlotSurface(L,varargin)
% Visualize the surface model
%
%   h = L.PlotSurface;
%   h = L.PlotSurface('opt1',val1,'opt2',val2,...);
%
% Available options:
%   - 'surface_idx', numbers. Give which surfaces to plot (1..L.nSurfaces).
%     Default is all surfaces.
%   - 'face_colors', color (see doc patch).
%     Defines face color per surface.
%   - 'edge_color', color | 'none'. Default is 'none'.
%   - 'face_alpha', numbers. Default is 0.5.
%   - 'noplot', false | true. Default is false
%     If true, then output is the patch data
%   - Other options given to patch command
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% default options
opt.surface_idx = 1:L.nSurfaces;
opt.face_colors = [0,0;1,0;0,1];
opt.edge_color = {'none'};
opt.face_alpha = [1,0.25];
opt.noplot = false;
others = {};

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else others = [others varargin(i:i+1)]; end
end

% adjust face colors
% if( size(opt.face_colors,2)~=numel(opt.surface_idx) )
%     opt.face_colors = repmat(opt.face_colors,1,numel(opt.surface_idx));
% end

% % adjust face alpha
% if( numel(opt.face_alpha)==1 )
%     opt.face_alpha = repmat(opt.face_alpha,1,numel(opt.surface_idx));
% end

% adjust edge color
if( numel(opt.edge_color)==1 )
    opt.edge_color = repmat(opt.edge_color,1,numel(opt.surface_idx));
end

if( ~opt.noplot ), hs = ishold; end

% for each surface
% apply transfor
FV.Vertices = L.surfacePoints;
for i=1:length(opt.surface_idx)
    si = opt.surface_idx(i);
    
    FV.Faces = L.GetSurfaceFaces(si);
    FV.FaceVertexCData = opt.face_colors(:,i);
    
    if( opt.noplot )
        h(i).FV = FV;
    else
    
        h(i) = patch(FV);
        set(h(i),others{:},'FaceColor',opt.face_colors(:,i),'EdgeColor',opt.edge_color{i},...
            'FaceAlpha',opt.face_alpha(i));
        hold on;
    end
    
end 

if( ~opt.noplot && ~hs ), hold off; end
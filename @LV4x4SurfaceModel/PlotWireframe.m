function h = PlotWireframe(L,varargin)
% Plot myocardial contour lines along longitudinal & circumferential
% directions along the nodal points.
%
%   h = L.PlotWireframe
%   h = L.PlotWireframe('opt1',val1,'opt2',val2,...)
%
% Available options:
%   - 'color', line color. Default is black.
%   - 'style', line style. Default is solid ('-').
%   - 'width', line width, Default is 2.0.
%   - 'long', indices for longitudinal lines. Set [] for no longitudinal lines. 
%     Default is 1:(L.nsamples(1)-1):L.nCirc-1 or [0 0.5*pi pi 1.5*pi]
%   - 'circ', indices for circular lines. Set [] for no circular line plots.
%     Default is 1:(L.nSamples(2)-1):L.nAzimuth
%   - 'with_nodes', true|false. Plot also the nodal points. Default is true.
%   - 'surfaces', string or cell of string. Default is {'epi','endo'}.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% default options
opt.color = 'k';
opt.style = '-';
opt.width = .5;
opt.circ = 1:(L.nSamples(2)-1):L.nAzimuth;
opt.circ = opt.circ(1:end-1); %skip mu=0
opt.long = 1:(L.nSamples(1)-1):L.nCirc-1; % thetas = [0 0.5pi pi 1.5pi]
opt.with_nodes = false;
opt.surfaces = {'epi','endo'};

% get option
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option %s',varargin{i}); end
end

ploto = {'Color','r','LineStyle',opt.style,'LineWidth',opt.width};
hs = ishold;

hold on;

% plot circular lines
P = L.GetCircumferentialPoints(opt.circ);
for i=1:size(P,4)
    h.circ(:,i) = plot3(P(:,:,1,i),P(:,:,2,i),P(:,:,3,i),ploto{:});
end

% plot longitudinal lines
P = L.GetLongitudinalPoints(opt.long);
for i=1:size(P,4)
    h.long(:,i) = plot3(P(:,:,1,i),P(:,:,2,i),P(:,:,3,i),ploto{:});
end

% plot the node points first
if( opt.with_nodes )
    cols = {'b','r'};
    for si=1:numel(opt.surfaces)
        N = L.GetNodalParameters(opt.surfaces{si});
        P = LV4x4SurfaceModel.Prolate2Cart(L.focalLength,N(:,3),N(:,2),N(:,1));
        
        % apply transformation
        P = [P ones(size(P,1),1)] * L.T;
        
        h.nodes(si) = plot3(P(:,1),P(:,2),P(:,3),'ks','MarkerFaceColor',cols{si});
    end

    % swap nodes
    h.nodes = fliplr(h.nodes);
end


if( ~hs ), hold off; end
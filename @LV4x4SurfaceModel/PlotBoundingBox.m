function h = PlotBoundingBox(L,varargin)
% Plot bounding box surrounding the LV model
%
%   h = L.PlotBoundingBox
%   h = L.PlotBoundingBox('opt1',val1,'opt2',val2,...)
%
% Available options:
%   - 'margin', number(s)
%     Set the margin off the outermost points on each dimension. Default is 5.
%
%   - 'color', line color. Default is black.
%   - 'width', line width. Default is 1.
%   - 'style', line style. Default is solid (see doc LineSpec).
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% default options
opt.margin = 5;
opt.color = 'k';
opt.width = 1;
opt.style = '-';

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ) opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option %s',varargin{i}); end
end

% adjust margin
if( numel(opt.margin)==1 ), opt.margin = repmat(opt.margin,1,3); end

% get bounding box
BB = L.GetBoundingBox('margin',opt.margin); 

% then plot
plotopts = {'Color',opt.color,'LineStyle',opt.style,'LineWidth',opt.width};
h(1) = line([BB(1:4,1); BB(1,1)],[BB(1:4,2); BB(1,2)],[BB(1:4,3); BB(1,3)],plotopts{:});
h(2) = line([BB(5:end,1); BB(5,1)],[BB(5:end,2); BB(5,2)],[BB(5:end,3); BB(5,3)],plotopts{:});
h(3) = line([BB(1,1); BB(5,1)],[BB(1,2); BB(5,2)],[BB(1,3); BB(5,3)],plotopts{:});
h(4) = line([BB(2,1); BB(6,1)],[BB(2,2); BB(6,2)],[BB(2,3); BB(6,3)],plotopts{:});
h(5) = line([BB(3,1); BB(7,1)],[BB(3,2); BB(7,2)],[BB(3,3); BB(7,3)],plotopts{:});
h(6) = line([BB(4,1); BB(8,1)],[BB(4,2); BB(8,2)],[BB(4,3); BB(8,3)],plotopts{:});

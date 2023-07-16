function BB = GetBoundingBox(L,varargin)
% Get bounding box around the LV model
%
%   BB = L.GetBoundingBox;
%   BB = L.GetBoundingBox('opt1',val1,'opt2',val2,...);
%
% Available options:
%   - 'margin', margin_xyz or [margin_x margin_y margin_z]
%     Set a number of pixels as the margin for each dimension. Default is [5 5 5].
%
% Author: Avan Suinesiaputra - Auckland Bionengineering Institute (2011)

if( isempty(L) ), error('Model is still empty.'); end

T = L.T;
L.T = eye(4);
L.T = L.T(:,1:3);

% default option
opt.margin = 5;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option %s',varargin{i}); end
end

% check margin
if( numel(opt.margin)==1 ), opt.margin = opt.margin * [1 1 1]; end

% get max & min of all points & adjust margin
minx = floor(min(L.surfacePoints)) - opt.margin;
maxx = ceil(max(L.surfacePoints)) + opt.margin;

% create BB
BB(1,:) = minx;
BB(2,:) = [minx(1) maxx(2) minx(3)];
BB(3,:) = [maxx(1) maxx(2) minx(3)];
BB(4,:) = [maxx(1) minx(2) minx(3)];
BB(5,:) = [minx(1) minx(2) maxx(3)];
BB(6,:) = [minx(1) maxx(2) maxx(3)];
BB(7,:) = maxx;
BB(8,:) = [maxx(1) minx(2) maxx(3)];

% put back L.T and apply T
BB = [BB ones(size(BB,1),1)] * T;
L.T = T;
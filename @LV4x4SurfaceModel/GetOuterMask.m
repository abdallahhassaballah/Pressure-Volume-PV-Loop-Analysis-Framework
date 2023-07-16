function [M,Lx] = GetOuterMask(L,varargin)
% Create an outer mask region, the output is a surface
%
%   [M,Lx] = L.GetOuterMask('opt1',val1,...);
%
% Available options:
%   - 'myo_scale'
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Insitute (2011)


M = [];
if( isempty(L) ), return; end


% default options
opt.myo_scale = 1.0;
opt.fl_scale = 1.0;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

% we assume that first part of lambdas is epi, the rest is endo
lambdas = L.lambdaParams;
lEpi = 1:(numel(lambdas)/2);
lEndo = (1+numel(lambdas)/2):numel(lambdas); 

% expand opt.scale times
meandist = mean(lambdas(lEpi) - lambdas(lEndo));
lambdas(lEpi) = lambdas(lEpi) + opt.myo_scale.*meandist;

% create new object and alter it's outer lambdas
Lx = LV4x4SurfaceModel('copyFrom',L);
Lx.lambdaParams = opt.myo_scale * lambdas;
Lx.focalLength = Lx.focalLength * opt.fl_scale;

% create the surface
M = Lx.GetClosedSurfaceTriangulation(Lx.nSurfaces);
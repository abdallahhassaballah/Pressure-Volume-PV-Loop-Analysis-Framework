function L = Create( focalLength, lambdaParams, muParams, thetaParams, varargin )
% Create an LV4x4SurfaceModel from scratch.
%
%   L = Create( focalLength, lambdaParams, muParams, thetaParams );
%   L = Create( focalLength, lambdaParams, muParams, thetaParams, ... );
%
% Inputs:  - focalLength = a scalar for the model's focal length,
%          - lambdaParams = 134 elements of global lambda parameters,
%          - muParams = 40 elements of global mu parameters,
%          - thetaParams = 40 elements of global theta parameters.
%
% Optional arguments:
%   - 'mapping', string or a struct.
%     Define GLM mapping to use. Default is 'CIM'.
%     Available mapping strings: 'CIM'.
%     If struct is used then it should have the following fields:
%       Lambda: [512x134] double
%           Mu: [128x40] double
%        Theta: [128x40] double
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% default options
opt.mapping = 'CIM';

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option.'); end 
end

% empty model
L = LV4x4SurfaceModel;

% mapping
if( isstruct(opt.mapping) )
    L.G2E = opt.mapping;
elseif( strcmpi(opt.mapping,'cim') )
    L = AssignCIMMapping(L);
else
    error('Unknown mapping %s',opt.mapping); 
end

% update parameters
L.SetParams(focalLength,lambdaParams(:),muParams(:),thetaParams(:));

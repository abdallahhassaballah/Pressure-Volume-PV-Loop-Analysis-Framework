function L = Create( focalLengths, lambdaParams, muParams, thetaParams, varargin )
% Create LVMM given the parameters
%
%   L = LVMM.Create( focalLengths, lambdaParams, muParams, thetaParams );
%   L = LVMM.Create( focalLengths, lambdaParams, muParams, thetaParams, ... );
%
% Available options:
%   - 'name', string
%   - 'ed', frame number
%   - 'es', frame number
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% check size
N = size(lambdaParams,2);
if( numel(focalLengths)==1 ), focalLengths = repmat(focalLengths,1,N); end

if( N ~= size(muParams,2) || N ~= size(thetaParams,2) || N ~= numel(focalLengths) )
    error('Invalid size.');
end

% default options
opt.name = '';
opt.ed = 0;
opt.es = 0;

% get options
for i=1:2:numel(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option %s',varargin{i}); end
end


% create
L = LVMM;
L.name = opt.name;

L.focalLengths = focalLengths;
L.lambdas = lambdaParams;
L.mus = muParams;
L.thetas = thetaParams;

L.ed = opt.ed;
L.es = opt.es;

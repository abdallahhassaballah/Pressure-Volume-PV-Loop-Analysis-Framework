function mass = GetMass(L,varargin)
% Calculate mass of the model (in gr)
%
%   mass = GetMass(L);
%   mass = GetMass(L,...);
%
% Optional arguments:
%   - 'coef', value. Default is 1.04 g/ml.
%
% Author: Avan Suinesiaputra - Centre Advanced Imaging, Univ. of Auckland (2012)

% default options
opt.coef = 1.04;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Option is not valid.'); end
end

mass = (L.GetEpiVolume - L.GetEndoVolume) * opt.coef;
function TF = isempty(L)
% Determine whether the model is empty or not
%
%   L.isempty or isempty(L)
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

TF = isempty(L.lambdaParams) || isempty(L.muParams) || ...
    isempty(L.thetaParams) || isempty(L.focalLength);
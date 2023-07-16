function [P,idx] = GetBasalPoints(L)
% Get the outer most basal points from an LV model
%
%   P = L.GetBasalPoints;
%   [P,idx] = L.GetBasalPoints;
%
% Input:   - L is an LV4x4SurfaceModel object
% Output:  - P are the 3D coordinate points of the basal points. Each
%            column contains points for each surface.
%          - idx is the same size of P indicating the element indices of
%            the corresponding basal points.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = []; idx = [];
if( isempty(L) ), return; end

% use the mapping, it's on the top
nSurfaceSamples = L.nAzimuth*L.nCirc;
for i=1:L.nSurfaces
    idx(i,:) = (i-1) * nSurfaceSamples + L.map(1,:);
    P(:,:,i) = L.surfacePoints(idx(i,:),:);
end


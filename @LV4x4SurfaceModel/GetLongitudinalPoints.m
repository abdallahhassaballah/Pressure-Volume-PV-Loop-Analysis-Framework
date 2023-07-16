function [P,idx] = GetLongitudinalPoints(L,ci)
% Get surface points that are on the same longitudinal angle.
%
%   P = L.GetLongitudinalPoints(ci);
%   [P,idx] = L.GetLongitudinalPoints(ci);
%
% Input:   - L is an LV4x4SurfaceModel object
%          - ci is the index of the circumferential component, which is in
%            the range of 1 to L.nCirc
% Output:  - P are the 3D coordinate points of the longitudinal points. Each
%            column contains points for each surface.
%          - idx is the same size of P indicating the element indices of
%            the corresponding longitudinal points.
%
% Note:
% - If numel(ci)==1, then size(P) = L.nAzimuth x 3 x L.nSurface
%   and size(idx) = L.nAzimuth x L.nSurface
% - If numel(ci)==N (N>1), then size(P) = L.nAzimuth x N x 3 x L.nSurface
%   and size(idx) = L.nAzimuth x N x L.nSurface
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = []; idx = [];
if( isempty(L) ), return; end

% get the correct elements for ci
nSS = L.nAzimuth * L.nCirc;
for i=1:L.nSurfaces
    idx(:,:,i) = (i-1) * nSS + L.map(:,ci);
    P(:,:,:,i) = reshape(L.surfacePoints(idx(:,:,i),:),[],numel(ci),3);
end

idx = squeeze(idx);
P = squeeze(P);
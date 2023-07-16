function [P,idx] = GetCircumferentialPoints(L,li)
% Get surface points that are on the same azimuthal position.
%
%   P = L.GetCircumferentialPoints(li);
%   [P,idx] = L.GetCircumferentialPoints(li);
%
% Input:   - L is an LV4x4SurfaceModel object
%          - li is the index of the azimuthal (longitudinal) component, which is in
%            the range of 1 to L.nAzimuth
% Output:  - P are the 3D coordinate points of the circumferential points. Each
%            column contains points for each surface.
%          - idx is the same size of P indicating the element indices of
%            the corresponding circumferential points.
%
% Note:
% - If numel(li)==1, then size(P) = L.nCirc x 3 x L.nSurface
%   and size(idx) = L.Circ x L.nSurface
% - If numel(li)==N (N>1), then size(P) = L.nCirc x N x 3 x L.nSurface
%   and size(idx) = L.nCirc x N x L.nSurface
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = []; idx = [];
if( isempty(L) ), return; end

% get the correct elements for li
nSS = L.nAzimuth * L.nCirc;
for i=1:L.nSurfaces
    idx(:,:,i) = (i-1) * nSS + L.map(li,:)';
    P(:,:,:,i) = reshape(L.surfacePoints(idx(:,:,i),:),[],numel(li),3);
end

%idx = squeeze(idx);
%P = squeeze(P);
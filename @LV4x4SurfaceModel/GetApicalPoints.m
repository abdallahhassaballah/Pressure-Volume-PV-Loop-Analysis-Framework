function [P,idx] = GetApicalPoints(L)
% Get the apical tip points of an LV model
%
%   P = L.GetApicalPoints;
%   [P, idx] = L.GetApicalPoints(i);
%
% Output: - P contains the 3D coordinates of the apical points
%         - idx are indices of the L.surfacePoints that are the apical
%           points. If multiple surfaces are supplied, idx becomes a
%           matrix.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = []; idx = [];
if( isempty(L) ), return; end

% get the apex (the bottom one)
idx = L.map(end,:);

% offset
nElmts = numel(L.map);

% multiply for L.nSurfaces times
idx = repmat(idx(:),1,L.nSurfaces) + (ones(numel(idx),1) * ((0:L.nSurfaces-1) * nElmts));

% find only the requested surface
for i=1:size(idx,2)
    P(:,:,i) = L.surfacePoints(idx(:,i),:);
end

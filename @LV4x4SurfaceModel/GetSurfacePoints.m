function [P,idx] = GetSurfacePoints(L,si)
% Get surface sample points at a particular surface index (1=endo .. L.nSurfaces=epi)
%
%   P = L.GetSurfacePoints(si)
%   [P,idx] = L.GetSurfacePoints(si)
%
% Input:  - si is surface index from 1 to L.nSurfaces or 'endo' or 'epi'
% Output: - P are the sampling points
%         - idx are the indices in the element matrix that correspond to P
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = []; idx = [];
if( isempty(L) ), return; end

if( ischar(si) )
    if( strcmpi(si,'endo') ), si = 1;
    elseif( strcmpi(si,'epi') ), si = L.nSurfaces;
    else error('Unknown surface %s',si); 
    end
else
    if( si<1 || si>L.nSurfaces ), error('Surface index is out of range.'); end
end


% get indices
idx = (si-1) * L.nAzimuth * L.nCirc + L.map;

% get P
P = L.surfacePoints(idx,:);
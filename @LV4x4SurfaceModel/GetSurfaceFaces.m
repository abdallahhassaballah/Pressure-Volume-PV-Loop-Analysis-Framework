function F = GetSurfaceFaces(L,si)
% Get the faces definition for a surface triangulation (see doc patch).
%
%   F = L.GetSurfaceFaces(si)
%
% Input:  - si is the surface index between 1 to L.nSurfaces,
%           or 'endo' or 'epi'
% Output: - F is Nx3 faces where N is the number of faces and each column
%           points to a L.surfacePoints row index that defines a vertex.
%           Each row is one triangle face.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

F = [];
if( isempty(L) ), return; end

if( ischar(si) )
    if( strcmpi(si,'endo') ), si = 1;
    elseif( strcmpi(si,'epi') ), si = L.nSurfaces;
    else error('Unknown surface %s',si); 
    end
else
    if( si<1 || si>L.nSurfaces ), error('Surface index is out of range.'); end
end

% If 2D mesh grid of the surface is defined by this indices:
%    1   2   3
%    4   5   6
%    7   8   9
%
% then the faces are define like this
% F1 = 1 4 2
% F2 = 2 5 3
% F3 = 4 7 5
% F4 = 5 8 6
% F5 = 2 4 5
% F6 = 3 5 6
% F7 = 5 7 8
% F8 = 6 8 9

elIds = [L.map L.map(:,1)]; % wrap around theta

% the top half
F1 = [ ...
    reshape(elIds(1:end-1,1:end-1)',[],1) ...
    reshape(elIds(2:end,1:end-1)',[],1) ...
    reshape(elIds(1:end-1,2:end)',[],1)];

% the bottom half
F2 = [ ...
    reshape(elIds(1:end-1,2:end)',[],1) ...
    reshape(elIds(2:end,1:end-1)',[],1) ...
    reshape(elIds(2:end,2:end)',[],1)];

% combine
F = [F1; F2];


% mind the surface index
F = (si-1) * numel(L.map) + F;
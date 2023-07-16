function [IN,ON] = InsideSurface(L,si,P,varargin)
% Determine whether points are inside a surface or not.
%
%   IN = L.InsideSurface(si,P,..);
%   [IN,ON] = L.InsideSurface(si,P,...);
%
% Input:  - si is the surface index between 1 to L.nSurfaces
%           or 'endo' or 'epi' string
%         - P is Nx3 point coordinates to test their location.
% Output: - IN is a vector of P indices that are inside a surface.
%         - ON is a vector of P indices that are on a surface.
%
% Available options:
%   - 'tol', number. Default is 1e-10.
%     Define surface tolerance to determine points on the surface.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

IN = []; ON = [];
if( isempty(L) ), return; end

if( ischar(si) )
    if( strcmpi(si,'endo') ), si = 1;
    elseif( strcmpi(si,'epi') ), si = L.nSurfaces;
    else error('Unknown surface %s',si); 
    end
else
    if( si<1 || si>L.nSurfaces ), error('Surface index is out of range.'); end
end

% create TriRep of the surface (closed)
FV0 = L.GetClosedSurfaceTriangulation(si);
T0 = TriRep(FV0.Faces,FV0.Vertices(:,1),FV0.Vertices(:,2),FV0.Vertices(:,3));

% centers & normal vectors
C0 = T0.incenters;
N0 = T0.faceNormals;

% find the closest point
np = size(P,1);
ns = size(C0,1);
[~,idxmin] = min(sqrt(((P(:,1)*ones(1,ns)) - (ones(np,1)*C0(:,1)')).^2 + ...
    ((P(:,2)*ones(1,ns)) - (ones(np,1)*C0(:,2)')).^2 + ...
    ((P(:,3)*ones(1,ns)) - (ones(np,1)*C0(:,3)')).^2),[],2);

% dot product these points with the normal vectors
dotV = dot(P-C0(idxmin,:),N0(idxmin,:)-C0(idxmin,:),2);

% INside is where dotV > 0
IN = find(dotV>0);

% ON the surface is where dotV == 0 (tolerance = 1e-6)
ON = find(abs(dotV)<1e-10);
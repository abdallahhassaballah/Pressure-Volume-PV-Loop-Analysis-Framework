function Ps = ProjectPointsToSurface(L,si,P)
% Project a set of points onto a surface
%
%   Ps = L.ProjectPointsToSurface(si,P);
%
% Input:  - si is the surface index from (1=endo) to (L.nSurfaces=epi).
%         - P are Nx3 points to be projected
% Output: - Ps are Nx3 projected points.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

Ps = [];
if( isempty(L) ), return; end

% check si
if( si<1 || si>L.nSurfaces )
    error('Invalid surface index.');
end

% get surface triangulation
Faces = L.GetSurfaceFaces(si);
Vertices = L.surfacePoints(Faces,:);
Faces = reshape(1:size(Vertices,1),size(Faces));

% create TriRep representation, normal faces and center of faces
Tri = TriRep(Faces,Vertices(:,1),Vertices(:,2),Vertices(:,3));
%NF = Tri.faceNormals;
CF = Tri.circumcenters;

% find the minimum distance from each P to the triangles
PCi = repmat(P,[1 1 numel(iNum)]) - permute(repmat(CF,[1 1 size(Pi,1)]),[3 2 1]);
[~,imin] = min(squeeze(sum(PCi.^2,2)),[],2);

